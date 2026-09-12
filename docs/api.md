# API Reference

For a runnable introduction, use [the quickstart](quickstart.md) or
[playback guide](playback.md). The top-level facade exports `Optimizer` and
`DenseArray`. Realized-array contracts and playback are explicit submodules.

## Optimization

::: dense_arrays.optimizer

## Dense-array results

::: dense_arrays.solution

## Sequence utilities

::: dense_arrays.sequence

## Realized arrays and playback

These contracts accept persisted sequence placements independently of the
optimizer. Producer adapters supply the translation from their own records.
See [playback authority](architecture/solution-playback.md) before interpreting
ordering or adding a renderer.

::: dense_arrays.realized

::: dense_arrays.playback
    options:
      members:
        - PlaybackPlan
        - PlaybackStep
        - PlaybackAuthority
        - OrderingStatus
        - ConstraintResult
        - PlaybackDocument
        - reconstruct_playback
        - dumps_realized_array
        - loads_realized_array
        - dumps_playback_plan
        - loads_playback_plan
        - render_playback_html
        - render_playback_collection_html
