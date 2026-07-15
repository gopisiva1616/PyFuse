# Conda recipe notes

This recipe is a starter template for conda-forge submission.

## Before submitting

1. Set metadata values:
- `home` URL in `meta.yaml`
- `recipe-maintainers` handles in `meta.yaml` (for conda-forge submission)

2. Build locally:

```bash
conda build conda-recipe
```

3. Smoke test package from local build.

4. Submit to conda-forge via staged-recipes or a feedstock.
