Format classes can now declare experimental metadata that they cannot read from the image
file, by listing PHIL paths in the ``missing_metadata`` class attribute or by calling
``declare_missing_metadata()`` from ``_start()``. The declared items are filled in from the
user's ``geometry.*`` overrides as the Format instance is constructed, and a
``MissingMetadataError`` is raised for anything left unsupplied. This removes the need for
such formats to invent dummy values that the user might silently fail to override.
