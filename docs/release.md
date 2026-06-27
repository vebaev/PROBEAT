# Release checklist

1. Confirm the worktree contains only reviewed release changes.
2. Generate and review Linux Conda lock files.
3. Verify curated-database provenance and licenses.
4. Run static, unit, smoke, and BV7/DPL5 regression tests.
5. Build the full Docker image and run `probeat doctor --profile full --strict`.
6. Run one real PE and one real ONT analysis through Streamlit.
7. Update `CHANGELOG.md`, `CITATION.cff`, and `manifests/release.yaml`.
8. Tag the reviewed commit and publish the immutable image digest in the
   GitHub Release.
