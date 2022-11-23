# Question about dissection variation
## Are the number of spots different across pos
## Are the percentage of each domain change across pos
## Per pos, is the cell composition of layers the same
## Cell composition change across samples



# # of Samples per pos ------------------------------------------------------------

fnl_dat |>
    group_by(pos) |>
    summarize(n_sample = n_distinct(Br_num),
              n_spots = n(),
              n_sp9_miss = sum(is.na(sp9)),
              n_sp16_miss = sum(is.na(sp16)))
## TODO: Are the number of spots different across pos
## Why missing dots
## TODO: with in each pos, is there particular sample that missing a lot
fnl_dat |>
    group_by(pos, Br_num) |>
    summarize(n_spots = n(),
              n_sp9_miss = sum(is.na(sp9)),
              n_sp16_miss = sum(is.na(sp16)),
              avg_sp9_miss = n_sp9_miss/n_spots,
              avg_sp16_miss = n_sp16_miss/n_spots)

## TODO: with in each pos, is there particular sample that missing a lot
fnl_dat |>
    group_by(pos, Br_num) |>
    summarize(n_spots = n(),
              n_sp9_miss = sum(is.na(sp9)),
              n_sp16_miss = sum(is.na(sp16)),
              avg_sp9_miss = n_sp9_miss/n_spots,
              avg_sp16_miss = n_sp16_miss/n_spots)


## Are the percentage of each domain change across pos

## Per pos, is the cell composition of layers the same

## Cell composition change across samples

fnl_dat |>
    group_by(pos, sample) |>
    summarize(n = n_distinct(Br_num))


# Sp9 Analysis ------------------------------------------------------------
fnl_dat |>
    group_by(pos, sp9) |>
    summarize(n_sample = n_distinct(Br_num),
              n_spots = n())
## TODO: Are the number of spots different across pos



# Sp16 Analysis ------------------------------------------------------------
