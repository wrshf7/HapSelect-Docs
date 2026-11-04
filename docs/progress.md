# Progress Reporting

Many HapSelect steps take minutes to hours on real datasets. HapSelect reports what it is doing and how far along it is in one of two ways:

- **Progress bars** (default) - one labelled bar per step, for people watching an R session or terminal.
- **Machine-readable events** - one JSON line per event on `stderr`, for other programs (pipelines, dashboards, web front ends) to follow a run. See [Machine-Readable Progress](#machine-readable-progress).

Progress reporting never changes any results.

---

## Progress Bars

<video class="method-video progress-video" src="../assets/progress-terminal.mp4" poster="../assets/progress-terminal-poster.jpg" autoplay loop muted playsinline></video>
<p class="method-video-caption progress-video">Progress bars from LD calculation through to simulation on the example data (sped up).</p>

Every long-running function reports as a named **step** with its own bar. The text beside a bar says what the step is currently doing (e.g., `Fold 3/5`, `Building blocks`, or the genetic algorithm's current fitness). Work that cannot be counted, such as waiting on BEAGLE to start, shows a spinner with the time taken. A step that stops with an error marks its bar as failed.

Bars are drawn when R is running interactively (e.g., RStudio) or in a terminal. When output is redirected to a file, as in most batch jobs, bars are switched off automatically and the same information is printed as plain text instead (for example, the genetic algorithm prints one line per iteration when `monitor = TRUE`).

!!! tip
    While a bar is being drawn, anything a step prints (such as the variance components printed for each cross-validation fold) is held back and printed once the step finishes, so it does not break up the bar. Nothing is lost.

### Turning progress bars off

```r
#for the rest of the session
options(progressr.enable = FALSE)
```

Or, before starting R (e.g., for a scheduled job), set the environment variable `R_PROGRESSR_ENABLE=false`.

### Other settings

| Setting | Default | Description |
|:---|:---|:---|
| `options(HapSelect.progress_interval = )` | `0.2` | Minimum number of seconds between progress updates. Must be set **before** HapSelect is loaded (e.g., in your `.Rprofile`) |
| `options(cli.progress_show_after = )` | `0` | Seconds before a bar is shown |
| `options(cli.progress_clear = )` | `FALSE` | If `TRUE`, finished bars are removed instead of left on screen |

### Using your own progress handler

HapSelect reports progress through the [progressr](https://progressr.futureverse.org/) package. If you have already set your own handlers with `progressr::handlers()` (e.g., in your `.Rprofile`) before HapSelect is loaded, HapSelect leaves them alone. You can also replace HapSelect's handler at any time:

```r
progressr::handlers("txtprogressbar")
```

Other handlers still show progress, but not the step names or status text, which only HapSelect's own handler knows about.

---

## Machine-Readable Progress

Setting the environment variable `HAPSELECT_MACHINE_LOGGING=1` **before** HapSelect is loaded replaces the progress bars with a stream of events on `stderr`, one JSON object per line:

```bash
HAPSELECT_MACHINE_LOGGING=1 Rscript my_analysis.R
```

```text
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:18.374Z","event":"step_started","step":"window_blocking","label":"Window Blocking"}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:18.385Z","event":"progress","step":"window_blocking","current":0,"total":21}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:18.390Z","event":"progress","step":"window_blocking","current":21,"total":21}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:18.391Z","event":"step_finished","step":"window_blocking","status":"succeeded"}
```

Nothing else about a run changes: results, printed output, warnings, and errors are the same as without it.

### Reading the stream

- Every event line starts with the prefix `@@HS ` (including the trailing space). Warnings, messages, and output from PLINK and BEAGLE share `stderr`, so keep only lines that start with the prefix and parse the rest of the line as JSON.
- Fields with no value are left out rather than sent as `null`.
- Numbers are written in full precision.
- Progress events are throttled to at most one every `HapSelect.progress_interval` seconds (0.2 by default), so not every increment is sent. Use `current` and `total` rather than counting events.

### Common fields

Every event has:

| Field | Type | Description |
|:---|:---|:---|
| `v` | string | Protocol version, which is the installed HapSelect version (e.g., `"2.0.0"`) |
| `at` | string | Time the event was sent, ISO 8601 in UTC with milliseconds (e.g., `"2026-11-02T04:14:18.374Z"`) |
| `event` | string | `"step_started"`, `"progress"`, `"step_finished"`, or `"metric"` |

### Events

| Event | Fields | Sent when |
|:---|:---|:---|
| `step_started` | `step`, `label` | A step begins. `label` is a human-readable name for display |
| `progress` | `step`, `current`, `total`, `message` | A step's count advances, or its status changes. All but `step` may be absent - see below |
| `step_finished` | `step`, `status`, `message` | A step ends. `status` is `"succeeded"` or `"failed"`, and a failed step carries the error in `message` |
| `metric` | `step`, `name`, `values` | A step produces a result worth charting live, such as one generation of the genetic algorithm. `values` is an object of named numbers |

**`progress` events**

- `current` and `total` are whole numbers, and `current` never exceeds `total`. `total` is absent when the amount of work is not known.
- `message` is the status text shown beside a progress bar (e.g., `"Fold 3/5"`).
- A status sent during work that is not counted (e.g., `"Running Beagle"`) repeats the step's last `current` and `total`, so it should not be read as the count starting over. Before a step has counted anything, a status has no counts at all.

```text
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.206Z","event":"progress","step":"beagle_phase","message":"Writing VCF"}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.953Z","event":"progress","step":"read_vcf","current":112,"total":112,"message":"Converting genotypes"}
```

**Nested steps**

Steps can run inside other steps. For example, `beagle_phase_geno()` reads BEAGLE's output with the VCF reader, so a `read_vcf` step starts and finishes inside the `beagle_phase` step. Track steps with a **stack**, not a single current step: `progress` and `metric` events belong to the innermost running step, named in their `step` field.

```text
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.199Z","event":"step_started","step":"beagle_phase","label":"Phasing (Beagle)"}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.947Z","event":"step_started","step":"read_vcf","label":"Reading VCF"}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.987Z","event":"step_finished","step":"read_vcf","status":"succeeded"}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.987Z","event":"step_finished","step":"beagle_phase","status":"succeeded"}
```

**Failed steps**

```text
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.988Z","event":"step_started","step":"graph_blocking","label":"Graph Blocking"}
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:28.988Z","event":"step_finished","step":"graph_blocking","status":"failed","message":"geno must be supplied for a graph_strategy: the method computes LD from the genotypes directly. Pass it to def_blocks(strategy, map, geno = ...)."}
```

### Steps

The `step` identifiers below are part of the protocol and are safe to key on. Labels are for display and may be reworded between versions.

| `step` | `label` | Reported by |
|:---|:---|:---|
| `read_vcf` | Reading VCF | `read_vcf_geno()`, `read_vcf_phased()`, and nested inside the `beagle_*_geno()` functions |
| `beagle_impute` | Imputation (Beagle) | `beagle_impute()`, `beagle_impute_geno()` |
| `beagle_phase` | Phasing (Beagle) | `beagle_phase_geno()` |
| `pairwise_ld` | Pairwise LD | `pairwise_ld()` |
| `plink_ld` | Pairwise LD (PLINK) | `plink_pairwise_ld()`, `plink_pairwise_ld_geno()` |
| `ld_blocking` | LD Blocking | `def_blocks()` with `ld_strategy()` |
| `window_blocking` | Window Blocking | `def_blocks()` with `window_strategy()` |
| `graph_blocking` | Graph Blocking | `def_blocks()` with `graph_strategy()` |
| `marker_effects` | Marker Effects | `solve_marker_effects()`, `create_marker_effects_file()` |
| `cross_validation` | Cross Validation | `cross_validation()` |
| `n_fold_cross_validation` | N-Fold Cross Validation | `n_fold_cross_validation()` |
| `local_gebv` | Local GEBV | `compute_local_GEBV()` |
| `haplotype_effects` | Haplotype Effects | `compute_haplotype_effects()` |
| `gebv_chunk_prep` | Preparing Blocks | Nested inside `local_gebv` / `haplotype_effects` when `parallel = TRUE` |
| `local_gebv_parent_selection` | Local GEBV Parent Selection | `local_gebv_parent_selection()` |
| `haplotype_parent_selection` | Haplotype Parent Selection | `haplotype_parent_selection()` |
| `local_gebv_simulation` | Local GEBV vs Truncation Simulation | `localGEBV_vs_TS_simulation()` |
| `haplotype_simulation` | Haplotype vs Truncation Simulation | `Haplotype_vs_TS_simulation()` |
| `block_var_test` | Haploblock Variance Test | Haploblock variance test |

### Metrics

| `name` | `values` | Sent |
|:---|:---|:---|
| `variance_components` | `marker_variance`, `additive_genetic_variance`, `residual_variance`, `h2` | Each time marker effects are solved, including once per cross-validation fold |
| `prediction_accuracy` | `accuracy` | By `create_marker_effects_file()`, after its `marker_effects` step has finished, so with no `step` field |
| `cv_fold` | `fold`, `accuracy` | Once per cross-validation fold |
| `cv_summary` | `folds`, `mean_accuracy` | At the end of cross-validation |
| `ga_generation` | `iter`, `best`, `mean`, `sd` | Once per generation of the genetic algorithm, with the fitness of the best solution and the mean and standard deviation across the population |

```text
@@HS {"v":"2.0.0","at":"2026-11-02T04:14:19.427Z","event":"metric","step":"marker_effects","name":"variance_components","values":{"marker_variance":3.67727338501975e-05,"additive_genetic_variance":0.0620438023570159,"residual_variance":0.165621257872503,"h2":0.272522284686402}}
```

### Example: following a run from Python

```python
import json
import os
import subprocess

env = dict(os.environ, HAPSELECT_MACHINE_LOGGING="1")
proc = subprocess.Popen(["Rscript", "my_analysis.R"], env=env,
                        stderr=subprocess.PIPE, text=True)

steps = []  #a stack, because steps can nest
for line in proc.stderr:
    if not line.startswith("@@HS "):
        continue  #a warning, message, or PLINK/BEAGLE output
    event = json.loads(line[len("@@HS "):])

    if event["event"] == "step_started":
        steps.append(event["step"])
        print("started:", event["label"])
    elif event["event"] == "step_finished":
        steps.remove(event["step"])
        print("finished:", event["step"], event["status"], event.get("message", ""))
    elif event["event"] == "progress" and "total" in event:
        print(f"  {event['step']}: {event['current']}/{event['total']}")
    elif event["event"] == "metric":
        print("  metric:", event["name"], event["values"])

proc.wait()
```

!!! tip
    With machine logging on, progress is reported even when R is not interactive, so it works under `Rscript`, job schedulers, and subprocesses. Setting `R_PROGRESSR_ENABLE=false` as well stops the `progress` events, but `step_started`, `step_finished`, and `metric` events are still sent.
