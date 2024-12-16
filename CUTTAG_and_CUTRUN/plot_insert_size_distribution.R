# args
args = commandArgs(T)
workDir = args[1]
runInfo = args[2]
libType = args[3]

setwd(workDir)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

if (libType == 'CUTTAG') {
  maxSize <- 1000
} else if (libType == 'CUTRUN') {
  maxSize <- 700
} else {
  message("Invalid library type!")
}

# all sample insert
insert_size <- list.files("./", pattern = '_insert_size.txt$', ignore.case = FALSE) %>%
  lapply(., function(x){read.table(x, header = FALSE, col.names = c(gsub("_insert_size.txt$","", x), "Size"), check.names = FALSE)}) %>%
  purrr::reduce(., dplyr::full_join, by = 'Size') %>% 
  dplyr::arrange(Size) %>%
  tidyr::pivot_longer(!Size, names_to = "Sample", values_to = "Freq")

size_ratio <- insert_size %>%
  group_by(Sample) %>%
  mutate(Ratio = Freq/sum(Freq, na.rm = TRUE))

# plot
n_sample <- unique(insert_size$Sample) %>% length()

p <- size_ratio %>% 
  ggplot(., aes(x = Size, y = Ratio)) +
  geom_vline(xintercept = c(180, 360, 540), lty = 2, col = "black", lwd = 0.3) +
  geom_line(aes(color = Sample), linewidth = 0.5, linetype = 1) +
  scale_x_continuous(limits = c(0, maxSize)) +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.title = element_blank()
  )

if (n_sample <= 12) {
  p <- p + theme(legend.position = "right")
} else {
  p <- p + theme(legend.position = "none")
}

# save
ggsave(filename = paste0(runInfo, "_all_sample_insert_size_distribution.pdf"), plot = p, width = 12, height = 5)

cat("Finish plot!\n")
