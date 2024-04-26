The template handling mechanism is extended so that a template with a
single ``#`` is expanded to match non-zero padded sequential numbers.
For example, ``image_#.cbf`` will match ``image_1.cbf``, ``image_2.cbf``,
..., ``image_10.cbf`` and so on.

Using a single ``#`` to match up to 10 images _within_ a zero-padded
sequence continues to work as before. For example,
``dials.import template=insulin_1_01#.img`` will match the files
``insulin_1_010.img``, ``insulin_1_011.img``, ..., ``insulin_1_019.img``,
and no others.
