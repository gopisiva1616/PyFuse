# Third-Party Notices

This project bundles third-party assets in the package config assets directory for HTML report rendering.

## Included assets and licenses

1. jQuery
- File: `src/pyfuse/config/jquery-3.6.0.min.js`
- Upstream: https://jquery.com/
- License: MIT
- Source confirmation: https://jquery.com/license/
- Evidence in bundled file header: `/*! jQuery v3.6.0 | ... | jquery.org/license */`

2. DataTables JavaScript
- File: `src/pyfuse/config/jquery.dataTables.min.js`
- Upstream: https://datatables.net/
- License: MIT (for DataTables 1.10+)
- Source confirmation: https://datatables.net/license/
- Evidence in bundled file header: `/*! DataTables 1.13.6 ... datatables.net/license */`

3. DataTables CSS
- File: `src/pyfuse/config/jquery.dataTables.min.css`
- Upstream: https://datatables.net/
- License: MIT (same DataTables distribution family)
- Source confirmation: https://datatables.net/license/

4. Inter font files
- Files:
  - `src/pyfuse/config/Inter-Regular.woff2`
  - `src/pyfuse/config/Inter-SemiBold.woff2`
- Upstream: https://github.com/rsms/inter
- License: SIL Open Font License 1.1
- Source confirmation:
  - https://rsms.me/inter/
  - https://raw.githubusercontent.com/rsms/inter/v4.1/LICENSE.txt

## Notes

- PyFuse remains licensed under GNU GPLv3 (see `LICENSE.md`).
- Third-party components retain their respective licenses.
- When redistributing PyFuse, keep this notice file and original asset headers intact.
