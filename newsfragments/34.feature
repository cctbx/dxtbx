Change dxtbx format registry to using entry points

dxtbx now discovers format classes during configuration time instead of
at runtime. Format classes can either be added into the dxtbx/format
directory as usual, registered by other python packages using the
'dxtbx.format' entry point, or installed by the user via the
'dxtbx.install_format' command.

To register format classes stored in ~/.dxtbx you need to run
'dxtbx.install_format -u' whenever you add or remove format classes.

Changes for library users:
* A number of registry lookup methods were deprecated or removed.
* Exceptions from format .understand() methods are no longer discarded.
  Similarly, when no matching format was found the datablock find_format()
  methods now return 'None' and no longer raise exceptions.
  In both cases the caller will need to deal with the situation appropriately.
* Format classes must be named 'Format*', and must inherit either from
  other format classes or from the top-level format class, 'Format'.
  Base classes must be given as their original name and must therefore not
  contain '.'s.
