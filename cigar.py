"""
cigar is a simple library for dealing with cigar strings. the most useful
feature now is soft-masking from left or right. This allows one to adjust
a SAM record only by changing the cigar string to soft-mask a number of bases
such that the rest of the SAM record (pos, tlen, etc.) remain valid, but
downstream tools will not consider the soft-masked bases in further analysis.


>>> c = Cigar('100M')
>>> len(c)
100
>>> str(c)
'100M'
>>> list(c.items())
[(100, 'M')]


>>> c = Cigar('20H20M20S')
>>> len(c)
40
>>> str(c)
'20H20M20S'
>>> list(c.items())
[(20, 'H'), (20, 'M'), (20, 'S')]

>>> c.mask_left(29).cigar, c.cigar
('20H9S11M20S', '20H20M20S')

>>> c = Cigar('10M20S10M')
>>> c.mask_left(10).cigar
'30S10M'
>>> c.mask_left(9).cigar
'9S1M20S10M'
>>> Cigar('10S').mask_left(10).cigar
'10S'
>>> Cigar('10H').mask_left(10).cigar
'10H'
>>> Cigar('10H').mask_left(11).cigar
'10H'
>>> Cigar('10H').mask_left(9).cigar
'10H'

>>> Cigar('1M10H').mask_left(9).cigar
'1S10H'

>>> Cigar('5M10H').mask_left(9).cigar
'5S10H'

>>> c = Cigar('1S1H1S5H1S5M10H')
>>> c.mask_left(9).cigar == c.cigar
True

>>> c = Cigar('1S1H1S5H1S5M10H')
>>> c.mask_right(9).cigar == c.cigar
True
>>> c.mask_right(11).cigar
'1S1H1S5H1S4M1S10H'


"""
from __future__ import print_function
from itertools import groupby
from operator import itemgetter

__version__ = "0.1.3"

class Cigar(object):
    read_consuming_ops = ("M", "I", "S", "=", "X")
    ref_consuming_ops = ("M", "D", "N", "=", "X")

    def __init__(self, cigar_string):
        self.cigar = cigar_string

    def items(self):
        if self.cigar == "*":
            yield (0, None)
            raise StopIteration
        cig_iter = groupby(self.cigar, lambda c: c.isdigit())
        for g, n in cig_iter:
            yield int("".join(n)), "".join(next(cig_iter)[1])

    def __str__(self):
        return self.cigar

    def __repr__(self):
        return "Cigar('%s')" % self

    def __len__(self):
        """
        sum of MIS=X ops shall equal the sequence length.
        """
        return sum(l for l, op,in self.items() \
                               if op in Cigar.read_consuming_ops)

    def reference_length(self):
        return sum(l for l, op in self.items() \
                               if op in Cigar.ref_consuming_ops)

    def mask_left(self, n_seq_bases, mask="S"):
        """
        Return a new cigar with cigar string where the first `n_seq_bases` are
        soft-masked unless they are already hard-masked.
        """
        cigs = list(self.items())
        new_cigs = []

        c, cum_len  = self.cigar, 0
        for i, (l, op) in enumerate(cigs):
            if op in Cigar.read_consuming_ops:
                cum_len += l
            if op == "H":
                cum_len += l
                new_cigs.append(cigs[i])
            elif cum_len < n_seq_bases:
                new_cigs.append(cigs[i])
            else:
                # the current cigar element is split by the masking.
                right_extra = cum_len - n_seq_bases
                new_cigs.append((l - right_extra, 'S'))
                if right_extra != 0:
                    new_cigs.append((right_extra, cigs[i][1]))
            if cum_len >= n_seq_bases: break
        else:
            pass

        new_cigs[:i] = [(l, op if op in "HS" else "S") for l, op in
                new_cigs[:i]]
        new_cigs.extend(cigs[i + 1:])
        return Cigar(Cigar.string_from_elements(new_cigs)).merge_like_ops()


    @classmethod
    def string_from_elements(self, elements):
        return "".join("%i%s" % (l, op) for l, op in elements if l !=0)


    def mask_right(self, n_seq_bases, mask="S"):
        """
        Return a new cigar with cigar string where the last `n_seq_bases` are
        soft-masked unless they are already hard-masked.
        """
        return Cigar(Cigar(self._reverse_cigar()).mask_left(n_seq_bases, mask)._reverse_cigar())


    def _reverse_cigar(self):
        return Cigar.string_from_elements(list(self.items())[::-1])

    def merge_like_ops(self):
        """
        >>> Cigar("1S20M").merge_like_ops()
        Cigar('1S20M')
        >>> Cigar("1S1S20M").merge_like_ops()
        Cigar('2S20M')
        >>> Cigar("1S1S1S20M").merge_like_ops()
        Cigar('3S20M')
        >>> Cigar("1S1S1S20M1S1S").merge_like_ops()
        Cigar('3S20M2S')
        """

        cigs = []
        for op, grps in groupby(self.items(), itemgetter(1)):
            cigs.append((sum(g[0] for g in grps), op))

        return Cigar(self.string_from_elements(cigs))


if __name__ == "__main__":
    import doctest
    doctest.testmod()
