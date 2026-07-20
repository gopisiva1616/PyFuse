# Resource Selection Behavior

PyFuse chooses resources in this order when running `annotator`:

1. `--resource_path` if provided
2. latest installed bundle for selected `--genome`

If no managed bundle exists, PyFuse exits and asks you to install resources.

## Managed Root Precedence

### Resources Subcommands (`install`, `list`, `path`, `verify`)

1. CLI `--resource-root`
2. `PYFUSE_RESOURCE_ROOT` environment variable
3. `resource_root` in settings YAML
4. platform default data directory

### Annotator Mode

1. `PYFUSE_RESOURCE_ROOT` environment variable
2. packaged settings key `resource_root`
3. platform default data directory

`annotator` does not expose `--resource-root`.

## Useful checks

```bash
pyfuse resources list
pyfuse resources path --genome grch37
pyfuse resources verify --genome grch37
```
