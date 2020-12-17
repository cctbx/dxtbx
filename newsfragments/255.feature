Allow format classes to be marked as ``@abstract``. This means that they
will be considered and returned by the Registry search if they are the
best match, but are intended to represent an incomplete "category" of
format class that other classes build on, so cannot be instantiated.
