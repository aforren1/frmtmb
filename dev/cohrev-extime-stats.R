## The paired example timings from dev/cohrev-extime.R, summarized.
## Round 1 only: rounds 2 and later in the same process reuse a warm
## tape and run at 0.2 s, which is not what R CMD check measures.
lane_ex <- c(6.41, 2.77, 1.94, 3.53, 6.57, 2.95, 3.05)
lane_cal <- c(0.910, 0.340, 0.200, 0.210, 0.940, 0.220, 0.270)
base_ex <- c(2.74, 2.33, 2.61, 6.75, 3.10, 3.15, 2.82)
base_cal <- c(0.390, 0.220, 0.250, 0.890, 0.950, 0.280, 0.320)
cat(sprintf("lane example CPU: min %.2f median %.2f max %.2f; over 5 s on %d of %d\n",
            min(lane_ex), median(lane_ex), max(lane_ex),
            sum(lane_ex > 5), length(lane_ex)))
cat(sprintf("base example CPU: min %.2f median %.2f max %.2f; over 5 s on %d of %d\n",
            min(base_ex), median(base_ex), max(base_ex),
            sum(base_ex > 5), length(base_ex)))
cat(sprintf("control block CPU, identical arithmetic in every process: %.2f to %.2f, a factor of %.1f\n",
            min(c(lane_cal, base_cal)), max(c(lane_cal, base_cal)),
            max(c(lane_cal, base_cal)) / min(c(lane_cal, base_cal))))
cat(sprintf("correlation of example CPU with the control, all 14 runs: %.3f\n",
            cor(c(lane_ex, base_ex), c(lane_cal, base_cal))))
cat(sprintf("runs with control above 0.85: example CPU %s\n",
            paste(c(lane_ex, base_ex)[c(lane_cal, base_cal) > 0.85],
                  collapse = " ")))
cat(sprintf("runs with control below 0.85: example CPU %s\n",
            paste(c(lane_ex, base_ex)[c(lane_cal, base_cal) < 0.85],
                  collapse = " ")))
cat(sprintf("\nWelch t on the arms, unnormalized: p = %.3f\n",
            t.test(lane_ex, base_ex)$p.value))
cat(sprintf("Welch t on example CPU per unit of control: p = %.3f\n",
            t.test(lane_ex / lane_cal, base_ex / base_cal)$p.value))
