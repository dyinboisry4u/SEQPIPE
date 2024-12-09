# args
args = commandArgs(T)
workDir = args[1]
runInfo = args[2]
libType = args[3]
umiLen = args[4]

setwd(workDir)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(dplyr))

umiLen <- as.numeric(umiLen)
if (libType == 'qPRO') {
  offset <- (umiLen+1)*2
} else if (libType == 'rPRO') {
  offset <- umiLen+1
} else {
  message("Invalid library type!")
}

# all sample insert
insert_size <- list.files("./", pattern = '.hist$', ignore.case = FALSE) %>%
  lapply(., function(x){read.table(x, header = FALSE, col.names = c("Size", gsub(".hist$","", x)), check.names = FALSE)}) %>%
  purrr::reduce(., dplyr::full_join, by = 'Size') %>% 
  dplyr::arrange(Size) %>%
  dplyr::mutate(Size = Size - offset) %>%
  tidyr::pivot_longer(!Size, names_to = "Sample", values_to = "Freq")

# plot
n_sample <- unique(insert_size$Sample) %>% length()
p <- insert_size %>% 
  dplyr::filter(Size >= 10) %>%
  ggplot(., aes(x = Size, y = Freq)) +
  geom_vline(xintercept = c(10, 20, 30, 40), lty = 2, col = "black", lwd = 0.3) +
  annotate("rect", xmin = 0, xmax = 20, ymin = -Inf, ymax = Inf, fill = '#9B7EBD', alpha = 0.4) +
  annotate("rect", xmin = 20, xmax = 30, ymin = -Inf, ymax = Inf, fill = '#D4BEE4', alpha = 0.2) +
  geom_point(aes(color = Sample), size = 0.8) +
  # scale_y_continuous(trans='log10') +
  scale_x_continuous(limits = c(0, 120), breaks = c(0, 10, 20, 30, 40, 60, 80, 100, 120)) +
  theme_light() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(colour = "black"),
    axis.ticks = element_line(colour = "black"),
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.title = element_blank()
  )
if (n_sample <= 24) {
  p <- p + theme(legend.position = "right")
} else {
  p <- p + theme(legend.position = "none")
}
# save
ggsave(filename = paste0(runInfo, "_all_sample_insert_size_distribution.pdf"), plot = p, width = 12, height = 8)

# degradation ratio
degrade_ratio <- insert_size %>% 
  dplyr::group_by(Sample) %>%
  dplyr::summarise(Rate = sum(Freq[Size>=10 & Size<=20], na.rm=TRUE)/sum(Freq[Size>=30 & Size<=40], na.rm=TRUE)) %>%
  dplyr::mutate(Rate = round(Rate, 4)) %>%
  tidyr::unite(Text ,c("Sample", "Rate"), sep = ': ')
# save
degrade_ratio %>% unlist() %>% readr::write_lines(., paste0(runInfo, "_all_sample_degradation_ratio.txt"))

cat("Finish plot!\n")
