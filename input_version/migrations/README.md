# Input file migrations

This directory contains machine-readable descriptions of changes to the 4C input file format,
used by [FourCIPP](https://github.com/4C-multiphysics/fourcipp)'s migration tool to
automatically upgrade old input files. Each file here corresponds to one `INPUT_VERSION` bump
and is named after the version it upgrades *to*, e.g. `1.1.0.yaml`.

For the full schema (available migration `type`s and their fields), see
[FourCIPP's migration README](https://github.com/4C-multiphysics/fourcipp/blob/main/src/fourcipp/migration/migrations/README.md).
Only the fields checked by `check-input-version-migration` (see
`utilities/four_c_python/src/four_c_ci/input_version_migration_schema.py`) are validated on
the 4C side; FourCIPP remains the source of truth for how migrations are actually applied.

## Contributing a migration

If your pull request changes, renames, removes, or otherwise breaks compatibility with an
existing input file parameter or section, add a migration fragment describing the change.
Because the target `INPUT_VERSION` is not known while the PR is open, use one of these
generic, reserved filenames instead of a concrete version number:

- `patch_upgrade.yaml` — small changes that only affect few users
- `minor_upgrade.yaml` — more severe changes that affect many users
- `major_upgrade.yaml` — complete change of the input file format (can probably not be captured by the migration tool)

Only one reserved fragment may exist at a time across all merged-but-unreleased PRs; if one
already exists, add your entries to it instead of creating a second one.

Once your PR is merged to `main`, an automated workflow renames the fragment to
`<new_version>.yaml` and updates `INPUT_VERSION` in a follow-up commit — you do not need to
do this yourself, and you do not need to know the resulting version number in advance.

### Example fragment (`minor_upgrade.yaml`)

```yaml
migrations:
  - id: add-new-optional-parameter
    type: parameter_added
    description: "Added the new optional parameter FOO to section BAR."
    path: ["BAR", "FOO"]
    value: 0.0
```
