# Same-snapshot derivative-tail replay (2026-11-18)

The comparison executable and all original inputs/results remain frozen at
`7321e6744932200176a0a6628b07e195d342071c`. These probes link the updated release
objects from the implementation worktree and read exactly the same saved old and
last-trial states. Object/source hashes, snapshot hashes and full diagnostic
summaries identify both versions. The new objects were built in a separate clean detached worktree at
`b0d33cb25ab2480b3eb0f5d37c18edc340b0f122`. The earlier intermediate replay
is retained in `../chexal_velocity_tail_replay_20261118`; this final replay
also includes the subsequently completed temporal-refinement failure.

| Chexal run | Before | After |
|---|---|---|
| SS3 N51, 42.1820094 s | All five Newton proposals fail point status 53 | All five have admissible cells and junction status 0 |
| OSC N51, 93.1000041 s | Two proposals status 71, three status 53 | Two remain 71; three become 0 |
| OSC N101, 94.2452856 s | One proposal status 71, four status 53 | One remains 71; four become 0 |
| OSC N51 dt=.0625, 77.3649936 s | All five Newton proposals fail point status 53 | All five have admissible cells and junction status 0 |
| SS1 N51, 52.7952655 s | Cells admissible; five MATCHED junction failures 1061 | Same distinction persists |

Here point status `53 = 30 + 10 + 13` identifies the closure derivative failure;
`71 = 10 + 61` is primitive recovery rejecting the inclined countercurrent state.
Junction `1061 = 1000 + 61` is the negative horizontal gas trace. No velocity,
root, spectrum or trace was clipped. The old stall near positive liquid velocity
`2.2586e-8 m/s` was numerical; it must not be labeled physical countercurrent flow.
Trials that actually leave the primary-source domain continue to fail.

The probe recomputes the linear solve with each object version. It compares the
same frozen states and five prescribed Newton factors, not bit-identical Newton
directions. This is a local diagnostic and does not establish long-run completion
or experimental validity. Fresh campaigns must use a separately frozen version.

Reproduction requires the saved run folders and release objects for each version:

```sh
python3 reproduce_snapshot.py --repo /path/to/frozen7321e67 --closure 5 --tag old_replay
python3 reproduce_snapshot.py --repo /path/to/frozen7321e67 --objects-repo /path/to/corrected-build --closure 5 --tag new_replay
python3 reproduce_snapshot.py --repo /path/to/frozen7321e67 --objects-repo /path/to/corrected-build --closure 5 --case D_OSC --label comparison_h_20261118 --tag new_replay
```

`before_after.json` supplies the concise mapping. Complete probe logs and summaries
are retained beside it; no source PDF content is included.
