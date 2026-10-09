---
title: Releases
description: Qualify immutable source tags and publish verifiable Python distributions.
---

# Release Dense Arrays

`pyproject.toml` owns the package version. Use immutable annotated `v<version>`
tags for reviewed releases: patches for compatible fixes, minor versions for
compatible additions, and major versions for incompatible public interfaces.
Persisted record schemas retain their own versions. Never move a published tag
or replace a published distribution.

The distribution is `dense-arrays`; imports remain `dense_arrays`. The historical
GitLab `paper` branch describes source installation of version 0.1.0. Preserve that frozen source and both authors' credits. Each downstream study
should identify the actual revision used for its results.

## Qualify and publish

1. Run the [development gate](../development.md#local-verification), inspect the
   wheel and source archive, and smoke-test a fresh installed wheel outside the
   checkout. Include a solved array and the rendering CLI; building alone is not
   an installation check.
2. Merge through a reviewed PR and require successful main-push CI for that exact
   commit. Check GitHub and PyPI for an unused version before creating its tag.
3. Configure a PyPI Trusted Publisher for `e-south/dense-arrays`, workflow
   `release.yaml`, environment `pypi`. Restrict the GitHub environment to `v*`
   tags. No API token belongs in source or workflow secrets.
4. Create the annotated version tag at the qualified main commit and publish its
   GitHub release. The release workflow verifies version, ancestry and exact CI,
   builds and tests without publishing credentials, then uploads from a separate
   OIDC publishing job.
5. Attach distributions and retained release evidence to the GitHub release.
   Evidence includes the source revision, successful CI run, runtime/build locks,
   artifact hashes and versioned citation. Workflow retention alone is temporary.
6. Compare published PyPI hashes with that evidence and install the exact version
   in a fresh environment before announcing availability. A failed upload is not
   a release; retry only the same qualified artifacts.

Downstream tools and papers pin a qualified version/artifact or immutable source
revision. Advancing that pin requires compatibility and numerical checks; the
package release itself does not qualify a paper's results.

See [PyPI Trusted Publishing](https://docs.pypi.org/trusted-publishers/creating-a-project-through-oidc/).

### Package-page images

The README is also the PyPI description. Keep its links absolute and use a
PNG banner at an immutable source commit or the matching `v<version>` tag.
Never use a relative asset path or a mutable branch for a release image. Keep
published tags and their assets; changing a new banner must not change older
release pages. When incrementing the version, update a version-bound README
image URL in the same change. The package tests enforce this relationship.

Before publishing, fetch the banner URL after the tag exists and compare its
bytes with the source PNG. Check the rendered PyPI page after upload. The PNG is
a package-page export; its editable SVG remains the artwork source.
