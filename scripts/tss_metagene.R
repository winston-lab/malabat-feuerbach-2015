library(psych)
library(tidyverse)
library(forcats)
library(ggforce)
library(ggthemes)

main = function(intable, trim_pct, downstream, ylabel, outpath){
    df = read_tsv(intable, col_names=c("group", "sample", "index", "position", "cpm")) %>% 
        mutate_at(vars(group), funs(fct_inorder(., ordered=TRUE)))
    nindices = max(df$index)
    df = df %>% group_by(group, sample, position) %>%
        summarise(trim_mean = winsor.mean(cpm, trim=trim_pct)) %>% 
        group_by(group, sample) %>%
        mutate_at(vars(trim_mean), funs(./max(trim_mean))) %>%
        ungroup()
    
    meta = ggplot(data = df, aes(x=position, y=trim_mean, group=sample, color=group)) +
            geom_vline(xintercept = 0, color="grey70") +
            geom_vline(xintercept = refsize/1000, color="grey70") +
            geom_hline(yintercept = 0, color="grey70") +
            geom_step(direction="vh", position=position_jitter(width=0), alpha=0.8) +
            scale_color_colorblind() +
            scale_x_continuous(breaks=c(0, refsize/1000), labels=c("TSS","CPS"),
                               expand=c(0,0)) +
            facet_zoom(xy= position > 0.2 & position < 1.85,
                       horizontal=FALSE, show.area=FALSE, zoom.size=1.5) +
            #facet_zoom(xy= position > -0.15 & position < .2,
            #           horizontal=FALSE, show.area=FALSE, zoom.size=1.5) +
            ylab("normalized counts") +
            ggtitle("mean sequencing coverage",
                    subtitle = paste(nindices, ylabel)) +
            theme_bw() +
            theme(text = element_text(size=12, face="bold", color="black"),
                  axis.text.y = element_text(size=10, color="black"),
                  axis.text.x = element_text(size=12, face="bold", color="black"),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(margin = margin(r=.5, unit="cm")),
                  legend.title = element_blank(),
                  legend.text = element_text(size=12),
                  plot.subtitle = element_text(face="plain"))
                  
    ggsave(outpath, meta, width=16, height=12, units="cm")
}

main(intable=snakemake@input[[1]],
     trim_pct=snakemake@params[["trim_pct"]],
     refsize=snakemake@params[["scaled_length"]],
     ylabel=snakemake@params[["ylabel"]],
     outpath=snakemake@output[[1]])
