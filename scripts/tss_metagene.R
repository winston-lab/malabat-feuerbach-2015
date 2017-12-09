library(psych)
library(tidyverse)
library(ggforce)
library(pals)

main = function(intable, trim_pct, refsize, outpath){
    df = read_tsv(intable,
                  col_names=c("group", "sample", "index", "position", "cpm")) %>%
            group_by(group, sample, position) %>%
            summarise(trim_mean = winsor.mean(cpm, trim=trim_pct)) %>% 
            group_by(group, sample) %>%
            mutate_at(vars(trim_mean), funs(./max(trim_mean))) %>%
            ungroup()
    
    meta = ggplot(data = df, aes(x=position, y=trim_mean, group=sample, color=group)) +
            geom_line() +
            scale_color_manual(values = brewer.set1(3)[2:3]) +
            facet_zoom(xy= position > 0.15 & position < 2,
                       horizontal=FALSE, show.area=FALSE, zoom.size=1.5) +
            scale_x_continuous(breaks=c(0, refsize/1000), labels=c("TSS","CPS")) +
            ylab("trimmed mean of normalized counts") +
            theme_bw() +
            theme(text = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_text(size=10),
                  axis.text.x = element_text(size=12, face="bold", color="black"),
                  legend.title = element_blank(),
                  axis.title.x = element_blank())
                  
    ggsave(outpath, meta, width=12, height=12, units="cm")
}

main(intable=snakemake@input[[1]],
     trim_pct=snakemake@params[["trim_pct"]],
     refsize=snakemake@params[["scaled_length"]],
     outpath=snakemake@output[[1]])
