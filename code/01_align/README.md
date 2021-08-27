# Create scripts

## Round 1

```R
sgejobs::job_single("round1", create_shell = TRUE, queue = "bluejay", memory = "35G", cores = 4L, task_num = 3, tc = 3, command = "module load cellranger/6.1.1")
```
