Cigar
=====

cigar is a simple library for dealing with cigar strings. the most useful
feature now is soft-masking from left or right. This allows one to adjust
a SAM record only by changing the cigar string to soft-mask a number of bases
such that the rest of the SAM record (pos, tlen, etc.) remain valid, but
downstream tools will not consider the soft-masked bases in further analysis.


```Python
>>> from cigar import Cigar

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

```

Installation
============

    pip install cigar
    conda install -c bioconda cigar
    
