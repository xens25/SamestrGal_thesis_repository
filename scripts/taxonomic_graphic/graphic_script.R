---
title: "Publication Figure 5"
output:
  html_notebook:
    code_folding: hide
    theme: united
    toc: yes
    toc_depth: 2
    toc_float: true
    toc_collapsed: false
theme: cosmo
---


# Libraries
Standard Import
```{r message=FALSE, warning=FALSE}
library(tidyverse)
library(cowplot)
library(scales)
library(ggpubr)
```

Special
```{r}
library(caret)
library(randomForest)
library(vip)
library(ggrepel)
library(grid)
library(tidytext)
library(tidymodels)
library(plotROC)
```

# Paths
```{r}
bin_dir = '../bin/'
data_dir = '../data/'
results_dir = '../results/'
```

# Custom Scripts / Theme
```{r}
source(paste0(bin_dir, 'theme_settings.R'))
source (paste0(bin_dir, 'utilities.R'))
source(paste0(bin_dir, 'species_competition_functions.R'))
source(paste0(bin_dir, 'distmat_utils.R'))
source(paste0(bin_dir, 'analysis_metadata.R'))
date <- format(Sys.time(), "%d%m%y")
```

# Import Tables
Metadata
```{r message=FALSE, warning=FALSE}
samples.metadata <- read_tsv(paste0(data_dir, 'samples.metadata.tsv'))
microbe.taxonomy <- read_tsv(paste0(data_dir, 'microbe.taxonomy.tsv')) 
microbe.metadata <- microbe.taxonomy %>% 
  right_join(., read_tsv(paste0(data_dir, 'microbe.metadata.tsv'))) %>% 
  mutate(gram.bool = ifelse(gram_stain == 'positive', T, 
                     ifelse(gram_stain == 'negative', F, NA)), 
         spores.bool = ifelse(spore_forming == 'spore-forming', T, 
                       ifelse(spore_forming == 'non-spore-forming', F, NA))) 

```

## Taxonomy
MetaPhlAn2 tables
combine with taxonomy
```{r message=FALSE, warning=FALSE}
mp_species.long <- microbe.taxonomy %>% 
  right_join(., read_tsv(paste0(data_dir, 'samples.mp_profiles.species.tsv')), 
             by = 'species') %>% 
  left_join(., samples.metadata, by = 'Name')
```

rCDI subset
annotate with metadata
```{r}
mp_species.rcdi.long <- 
  mp_species.long %>%
  filter(Study %in% c(rcdi.studies)) 
```

## Case Summary
SameStr Case-Summary table
```{r message=FALSE, warning=FALSE}
sstr_cases <- read_tsv(paste0(data_dir, 'samples.case_summary.tsv')) %>% 
  left_join(., microbe.taxonomy, by = 'species')
```

# Set Figure
```{r}
fig = paste0('Fig_5.')
```

# Strain-level Source (rCDI)
Strain-level plot of per-Case post-FMT rel. abund. by source

## Format Data
```{r fig.width=5, fig.height=5, message=FALSE, warning=FALSE, paged.print=FALSE}
sstr_cases.rcdi.metrics <-
  sstr_cases %>% 
  
  filter(Study %in% c('ALM', 'FRICKE')) %>% 
  filter(Case_Name %in% cases.full) %>% 

  tag_species_competition(.) %>% 

  mutate(n = 1) %>% 
  mutate(source = ifelse(grepl(species, pattern = 'unclassified'),  'Unclassified', source)) %>% 
  
  mutate(source = case_when(
  analysis_level == 'species' & source == 'self' ~ 'Self Sp.',
  analysis_level == 'species' & source == 'donor' ~ 'Donor Sp.',
  analysis_level == 'species' & source == 'unique' ~ 'Unique Sp.',
  T ~ source
  )) %>% 
  
  group_by(Study_Type, Case_Name, source, Days_Since_FMT.post, fmt_success.label) %>% 
  summarize(n = sum(n, na.rm = T),
            rel_abund = sum(rel_abund.post, na.rm = T) / 100) %>% 
  group_by(Case_Name, Days_Since_FMT.post) %>% 
  mutate(f = n / sum(n, na.rm = T)) %>%
  ungroup() %>% 

  group_by(Case_Name) %>% 
  mutate(source = case_when(source == 'Unclassified' ~ 'Unclassified Sp.',
                  source == 'same_species' ~ 'Same Sp.', 
                  source == 'unique' ~ 'Unique to Time Point',
                  source == 'self' ~ 'Recipient / Initial Sample',
                  source == 'donor' ~ 'Donor',
                  source == 'both' ~ 'Coexistence', 
                  T ~ source)) %>% 
  
  pivot_wider(names_from = 'source', 
               values_from = c('rel_abund','f', 'n'), names_sep = '___') %>% 
  mutate_at(.vars = vars(contains('___')),
            .funs = funs(replace_na(., 0))) %>%
  pivot_longer(cols = contains("___"), 
               names_to = c("metric", "source"),
               names_sep = '___', values_drop_na = F)
```

## Bar Chart, per Case
```{r fig.width=10}
pseudo = 0.5
bb = c(pseudo, 2, 7, 14, 35, 84, 365)

strain_order <- c('Donor', 'Donor Sp.', 'Coexistence', 'Recipient / Initial Sample', 'Self Sp.', 'Unique to Time Point', 'Unique Sp.', 'Same Sp.')
strain_colors <- c(colors.discrete[c(6,1,4, 8,3, 7,2, 5,10)])
use_vars = list(strain_order, strain_colors)

plot.rcdi <- 
  sstr_cases.rcdi.metrics %>% 
  filter(metric == 'rel_abund') %>% 
  ungroup() %>% 
  mutate(ordering = -as.numeric(as.factor(Case_Name)) + (Days_Since_FMT.post/1000)) %>% 
  ggplot(aes(
    fct_reorder(paste0('D', Days_Since_FMT.post, ' | ', str_remove(Case_Name, pattern = 'Case_')), ordering),
    y = value, 
    fill = fct_relevel(source, use_vars[[1]]))) + 

  geom_bar(stat = 'identity', position = position_fill(), width = 1, col = 'black') +
  theme_cowplot() + 
  scale_y_continuous(labels = percent_format(),
                     breaks = c(1, .75, .5, .25, 0),
                     expand = c(0,0)) +
  scale_fill_manual(values = strain_colors, 
                    
                    guide = guide_legend(reverse = TRUE)) + 
  labs(y = 'Rel. Abund.') + 
  theme(
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 10),
        legend.title=element_blank(),
        legend.position = 'top') 
    guides(fill=guide_legend(nrow = 2))

plot.rcdi
```

Exporting plots
```{r}
output_name = 'PostSource.Cases.rCDI.Strain'
ggsave(plot.rcdi + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 9, height = 3.5,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
```

## Area Chart, Across Cases
rCDI - Area chart showing post-FMT rel. abund. of sets of co-occurring strains in recipients, donors, and post-FMT patients based on glm 2nd binomial model across successful rCDI cases.
We only have sparse data and single cases for later time points, limiting to <= 84 days
```{r fig.width=3, fig.height=5, message=FALSE, warning=FALSE, paged.print=FALSE}
bb = c(2, 7, 14, 35, 84, 365)

strain_order <- c('Donor', 'Donor Sp.', 'Coexistence', 'Recipient / Initial Sample', 'Self Sp.', 'Unique to Time Point', 'Unique Sp.', 'Same Sp.')
strain_colors <- c(colors.discrete[c(6,1,4, 8,3, 7,2, 5,10)])
use_vars = list(strain_order, strain_colors)

plot <-

  sstr_cases.rcdi.metrics %>%   
  
  filter(fmt_success.label == 'Resolved', Days_Since_FMT.post <= 84) %>% 
  filter(metric != 'n') %>% 
  ggplot(aes(x = Days_Since_FMT.post,  
             y = value, 
             fill = fct_relevel(source, use_vars[[1]]))
         ) +  
  
  stat_smooth(geom = 'area',
              method = 'glm',  colour = 'black', size = 0.25,
              position = 'fill',
              # formula = y ~ poly(x, 2),
              method.args=list(family='binomial')
              ) +
  
   facet_grid(scales = 'free_x', space = 'free_x',
   cols = vars(''),
   rows = vars(ifelse(grepl(metric, pattern = '^f'), 
                                'Spp. Fraction', 
                                'Rel. Abund.'))) +
  theme_cowplot() + 
  scale_y_continuous(labels = percent_format(),
                     breaks = c(1, .75, .5, .25, 0),
                     expand = c(0,0), 
                     ) +
  scale_x_continuous(trans = pseudo_log_trans(),
                     breaks = bb,
                     expand = c(0,0)) + 
  scale_fill_manual(values = use_vars[[2]], 
                    guide = guide_legend(reverse = TRUE)) + 
  labs(x = 'Days Since FMT\n') + 
  theme(axis.title.y = element_blank(), 
        plot.margin = unit(c(0, 0, 0, 0), "cm"), 
        strip.background = element_blank(), 
        strip.text = element_text(size = 14), 
        panel.spacing.y = unit(.5, "cm"),
        # panel.spacing.x = unit(.35, "cm"),
        legend.title=element_blank(),
        legend.position = 'top', 
        legend.margin=ggplot2::margin(l = 1.5, unit='cm')) + 
    guides(fill=guide_legend(ncol = 2))
plot + theme(legend.position = 'none')
```

Histogram of available samples for modelling
```{r fig.width=5, fig.height=2, message=FALSE, warning=FALSE, paged.print=FALSE}
plot.histogram <-
  sstr_cases.rcdi.metrics %>%   
  
  filter(fmt_success.label == 'Resolved', Days_Since_FMT.post <= 84) %>% 
  filter(metric == 'f') %>% 
  
  group_by(Study_Type, Case_Name, Days_Since_FMT.post) %>% 
  summarize(value = sum(value, na.rm = T)) %>% 
  ggplot(aes(Days_Since_FMT.post + pseudo)) + 
  geom_histogram(fill = 'black', bins = 10) +
  facet_grid(scales = 'free_x', space = 'free_x', 
   cols = vars(fct_relevel(Study_Type, 'rCDI', 'Control')),
   rows = vars('n')) +
  scale_x_continuous(trans = pseudo_log_trans(0.5, 10), 
                     breaks = c(2, 7, 14, 35, 84, 365),
                     expand = c(0,0)) + 
  scale_y_continuous(#trans = 'log10', 
                     breaks = c(1, 5, 10, 100), 
                     expand = c(0,0),
                     ) +
  theme_cowplot() + 
  theme(strip.background = element_blank(), 
        strip.text.x = element_blank(), 
        axis.title = element_blank(), 
        axis.line.x = element_line(size = 1),
        axis.line.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title=element_blank())
plot.histogram
legend = cowplot::get_legend(plot)
grid.newpage()
grid.draw(legend)
```

Exporting plot
```{r}
output_name = 'PostSource.AcrossCases.rCDI.Strain'

ggsave(plot + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 3, height = 5,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
ggsave(plot.histogram + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 3, height = .6,
       filename = paste0(results_dir, fig, output_name, '.Histogram','.Plot.pdf'))
ggsave(legend, 
       device = 'pdf', dpi = 300, width = 5, height = 2,
       filename = paste0(results_dir, fig, output_name, '.Legend.pdf'))
```


Persistence / Engraftment within time periods:
```{r}
sstr_cases.rcdi.metrics %>% 
  filter(fmt_success.label == 'Resolved', metric == 'rel_abund') %>% 
  filter(Days_Since_FMT.post > 0, Days_Since_FMT.post <= 7) %>% 
  filter(source %in% c('Donor', 'Recipient / Initial Sample')) %>% 
  group_by(source) %>% 
  summarize(mean = mean(value) * 100, 
            sd = sd(value) * 100, 
            .groups = 'drop')

sstr_cases.rcdi.metrics %>% 
  filter(fmt_success.label == 'Resolved', metric == 'rel_abund') %>% 
  filter(Days_Since_FMT.post > 70, Days_Since_FMT.post <= 84) %>% 
  filter(source %in% c('Donor', 'Recipient / Initial Sample')) %>% 
  group_by(source) %>% 
  summarize(mean = mean(value) * 100, 
            sd = sd(value) * 100, 
            .groups = 'drop')
```
Donor -> 42.5 ± 30.3  to  26.5 ± 21.9
Recipient -> 18.9 ± 22.3  to  4.9 ± 9.0


## Grouped by Case
```{r fig.width=25, fig.height=5}
sstr_cases.rcdi.metrics.format <-
  sstr_cases.rcdi.metrics %>% 
  
  mutate(fmt_success.tag = ifelse(fmt_success.label == 'Resolved', 'S', 'F')) %>% 
  filter(metric == 'rel_abund') %>% 
  
  mutate(ALM = convert_to_alm_cases(Case_Name)) %>% 
  mutate(Study = str_split_fixed(Case_Name, pattern = '_', 2)[,1]) %>% 
  mutate(Case_Number = ifelse(Study == 'ALM', ALM, str_split_fixed(Case_Name, pattern = 'Case_', 2)[,2]))

sstr_cases.rcdi.metrics.format %>% 
  ggplot(aes(x = fct_reorder(paste0(Days_Since_FMT.post), Days_Since_FMT.post, .desc = F), 
             y = value, 
             fill = fct_relevel(source, use_vars[[1]]))
         ) +  
  geom_bar(stat = 'identity', show.legend = F, width = 1, col = 'black', position = position_fill()) + 
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(labels = percent_format_signif, expand = c(0,0)) + 
  scale_fill_manual(values = strain_colors) + 
  theme_cowplot() + 
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        
        legend.title = element_blank()) + 
  facet_grid(cols = vars(fmt_success.label, Study, Case_Number),
                         scales = 'free', space = 'free')

```

```{r}
sstr_cases.rcdi.metrics.out <- 
  sstr_cases.rcdi.metrics.format %>% 
  mutate(source = ifelse(source == 'Recipient / Initial Sample', 'Self', 
                         str_split_fixed(source, pattern = ' ', 2)[,1])) %>% 
  select(-metric, -ALM, -Study_Type, -fmt_success.tag) %>% 
  ungroup() %>% 
  pivot_wider(names_from = 'source', values_from = 'value', values_fn = list(value = sum)) %>% 
  mutate(Total = Donor + Same + Self + Unclassified + Unique + Coexistence, 
         DonorRecipientRatio = log10((Donor + 1e-6) / (Self + 1e-6))) 

write_tsv(sstr_cases.rcdi.metrics.out, paste0(results_dir, fig, 'FMT_Events.tsv'))
```

# Shared Strains per Taxa

## Format
Summarize per genus
```{r}
sstr_cases.study <- 
  sstr_cases %>%
  filter(Days_Since_FMT.post <= 373) %>% # same time frame rCDI / control
  filter(Study_Type == 'Control' | fmt_success) %>% 
  
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 

  mutate(source = case_when(
    analysis_level == 'species' & source == 'same_species' & Study_Type == 'Control' ~ 'self',
    analysis_level == 'species' & source == 'self' ~ 'Recipient / Initial Sample Sp.',
    analysis_level == 'species' & source == 'donor' ~ 'Donor Sp.',
    analysis_level == 'species' & source == 'unique' ~ 'Unique to Time Point Sp.',
    T ~ source
  )) %>% 
  
  mutate(source = case_when(source == 'Unclassified' ~ 'Unclassified Sp.',
                  source == 'same_species' ~ 'Same Sp.', 
                  source == 'unique' ~ 'Unique to Time Point',
                  source == 'self' ~ 'Recipient / Initial Sample',
                  source == 'donor' ~ 'Donor',
                  source == 'both' ~ 'Coexistence',
                  T ~ source)) %>% 
  
  filter(!source %in% c('Same Sp.') | Study_Type == 'Control') %>%
  
  mutate(genus = case_when(
    genus == 'Clostridiaceae_noname' ~ 'Clostridiales_noname',
    T ~ genus))

sstr_cases.species <- 
  sstr_cases.study %>%
  
  mutate(n = 1) %>% 
  group_by(Study_Type, source, species) %>% 
  summarize(n = sum(n, na.rm = T), .groups = 'drop') %>% 
  
  pivot_wider(names_from = 'Study_Type', values_from = 'n') %>%  
  mutate_at(.vars = vars(rCDI, Control),
            .funs = funs(replace_na(., 0))) %>% 
  
  mutate(source.simple = gsub(source, pattern = ' Sp.', replacement = ''))

sstr_cases.genus <-
  sstr_cases.study %>% 
  
  mutate(n = 1) %>% 
  group_by(Study_Type, source, genus) %>% 
  summarize(n = sum(n, na.rm = T), .groups = 'drop') %>% 
  
  pivot_wider(names_from = 'Study_Type', values_from = 'n') %>%  
  mutate_at(.vars = vars(rCDI, Control),
            .funs = funs(replace_na(., 0))) %>% 
  
  mutate(source.simple = gsub(source, pattern = ' Sp.', replacement = ''))
```

## Feature rates
Fraction of Taxa from shown post-FMT samples in which species have a feature
```{r}
genus_feature_rates <- 
  mp_species.rcdi.long %>%
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 
  filter(Sample_Type %in% c('post'), fmt_success) %>% 
  left_join(., microbe.metadata) %>%
  group_by(genus) %>% 
  summarize_at(.vars = vars(habit.oral, oxygen.tolerant, gram.bool, spores.bool), 
               .funs = funs(replace_na(
                 sum(. == T, na.rm = T) / sum(!is.na(.)), 
                                       0)))

species_feature_rates <- 
  mp_species.rcdi.long %>%
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 
  filter(Sample_Type %in% c('post')) %>% 
  left_join(., microbe.metadata) %>%
  group_by(species) %>% 
  summarize_at(.vars = vars(habit.oral, oxygen.tolerant, gram.bool, spores.bool), 
               .funs = funs(replace_na(
                 sum(. == T, na.rm = T) / sum(!is.na(.)), 
                                       0)))
```

## Taxa importance
Measured by observed events in each genus
```{r}
genus_importance <-
  sstr_cases.genus %>% 
  group_by(genus) %>% 
  summarize_at(.vars = vars(rCDI, Control), 
               .funs = funs(sum(., na.rm = T))) %>% 
  ungroup() %>% 
  mutate(total = rCDI + Control) %>% 
  filter(rCDI >= 5, Control >= 3) %>% 
  arrange(desc(rCDI))

species_importance <-
  sstr_cases.species %>% 
  group_by(species) %>% 
  summarize_at(.vars = vars(rCDI, Control), 
               .funs = funs(sum(., na.rm = T))) %>% 
  ungroup() %>% 
  mutate(total = rCDI + Control) %>% 
  filter(rCDI >= 5, Control >= 3) %>% 
  arrange(desc(rCDI))
```

## Transmission Rates
Control persistence and rCDI Engraftment rates for genera and species
```{r fig.width=8, fig.height=15, message=FALSE, warning=FALSE, paged.print=FALSE}
pseudo = 1e-10
genus_transmission_rates <- 
  sstr_cases.genus %>% 
  
  filter(genus %in% genus_importance$genus) %>% 
  
  pivot_longer(names_to = 'Study_Type', values_to = 'n', cols = c('Control','rCDI')) %>% 

  # summarize engraftment (strain + species)
  group_by(Study_Type, source.simple, genus) %>% 
  summarize(n = sum(n, na.rm = T)) %>% 
  
  group_by(Study_Type, genus) %>% 
  mutate(f = n / sum(n, na.rm = T)) %>% 
  
  # filter control-persistence and rcdi-engraftment
  mutate(label = case_when(
    source.simple == 'Recipient / Initial Sample' & Study_Type == 'Control' ~ 'Control.Persistence', 
    source.simple == 'Recipient / Initial Sample' & Study_Type == 'rCDI' ~ 'rCDI.Persistence', 
    source.simple == 'Donor' & Study_Type == 'rCDI' ~ 'rCDI.Engraftment',
    source.simple == 'Unique to Time Point' & Study_Type == 'Control' ~ 'Control.New')) %>% 
  filter(!is.na(label)) %>% 
  ungroup() %>% 
  
  # pivot
  select(genus, f, n, label) %>% 
  pivot_wider(names_from = 'label', 
              values_from = c('n','f'), 
              values_fill = list(n = 0, f = 0),
              names_sep = '.') %>% 
  
  # calculate percent ranks (pseudo), rescale to 0-1
  mutate(Control.Persistence.rank = rescale(dense_rank(f.Control.Persistence), to = c(pseudo, 1 - pseudo)), 
         rCDI.Engraftment.rank = rescale(dense_rank(f.rCDI.Engraftment), to = c(pseudo, 1 - pseudo)),
         rCDI.Persistence.rank = rescale(dense_rank(f.rCDI.Persistence), to = c(pseudo, 1 - pseudo)),
         Control.New.rank = rescale(dense_rank(f.Control.New), to = c(pseudo, 1 - pseudo)))

species_transmission_rates <- 
  sstr_cases.species %>% 
  
  filter(species %in% species_importance$species) %>% 
  
  pivot_longer(names_to = 'Study_Type', values_to = 'n', cols = c('Control','rCDI')) %>% 

  # summarize engraftment (strain + species)
  group_by(Study_Type, source.simple, species) %>% 
  summarize(n = sum(n, na.rm = T)) %>% 
  
  group_by(Study_Type, species) %>% 
  mutate(f = n / sum(n, na.rm = T)) %>% 
  
  # filter only control-persistence and rcdi-engraftment
  mutate(label = case_when(
    source.simple == 'Recipient / Initial Sample' & Study_Type == 'Control' ~ 'Control.Persistence', 
    source.simple == 'Recipient / Initial Sample' & Study_Type == 'rCDI' ~ 'rCDI.Persistence', 
    source.simple == 'Donor' & Study_Type == 'rCDI' ~ 'rCDI.Engraftment')) %>% 
  filter(!is.na(label)) %>% 
  ungroup() %>% 
  
  # pivot
  select(species, f, n, label) %>% 
  pivot_wider(names_from = 'label', values_from = c('n','f'), names_sep = '.') %>% 
  mutate_at(.vars = vars(everything(), -species), 
            .funs = funs(replace_na(., 0))) %>% 
  
  # calculate percent ranks (pseudo), rescale to 0-1
  mutate(Control.Persistence.rank = rescale(dense_rank(f.Control.Persistence), to = c(pseudo, 1 - pseudo)), 
         rCDI.Engraftment.rank = rescale(dense_rank(f.rCDI.Engraftment), to = c(pseudo, 1 - pseudo)))
```

```{r}
output_name = 'EngraftmentPersistence.Species'
species_transmission_rates %>% 
  rename_at(.vars = vars(starts_with('n.')), 
            .funs = funs(paste0(gsub(., pattern = '^n.', replacement = ''), '.count'))) %>% 
  rename_at(.vars = vars(starts_with('f.')), 
            .funs = funs(paste0(gsub(., pattern = '^f.', replacement = ''), '.fraction'))) %>% 
  select(species, starts_with('Control'), starts_with('rCDI')) %>% 
  write_tsv(., paste0(results_dir, fig, output_name, '.Table.tsv'))

output_name = 'EngraftmentPersistence.Genus'
genus_transmission_rates %>% 
  rename_at(.vars = vars(starts_with('n.')), 
            .funs = funs(paste0(gsub(., pattern = '^n.', replacement = ''), '.count'))) %>% 
  rename_at(.vars = vars(starts_with('f.')), 
            .funs = funs(paste0(gsub(., pattern = '^f.', replacement = ''), '.fraction'))) %>% 
  select(genus, starts_with('Control'), starts_with('rCDI')) %>% 
  write_tsv(., paste0(results_dir, fig, output_name, '.Table.tsv'))
```


## Per Genus
### Format
```{r fig.width=8, fig.height=15, message=FALSE, warning=FALSE, paged.print=FALSE}
sstr_cases.genus.fraction <-
  
  sstr_cases.genus %>% 
  
  ## filter genera
  filter(genus %in% genus_importance$genus) %>% 
  pivot_longer(names_to = 'Study_Type', values_to = 'n', cols = c('Control','rCDI')) %>% 

  ## calculate fraction of events
  group_by(Study_Type, genus, source.simple) %>%
  mutate(source.n = sum(n, na.rm = T)) %>% 
  
  group_by(Study_Type, genus) %>%
  mutate(total = sum(n, na.rm = T), 
         f = n / total,
         source.f = source.n / total) %>% 
  
  ## get genus-level donor-derived for ordering
  group_by(genus) %>% 
  mutate(donor_derived.genus = max(ifelse(Study_Type == 'rCDI' & 
                                          source.simple == 'Donor', 
                                          source.f, 0), na.rm = T)) %>% 
  ungroup()  


sstr_cases.genus.fraction.annotated <- 
  sstr_cases.genus.fraction %>% 
  
  # add genus metrics for present species
  left_join(., genus_feature_rates) %>%
  left_join(., genus_transmission_rates) %>%
  
  pivot_longer(names_to = 'feature', values_to = 'value', 
               cols = c(gram.bool, habit.oral, 
                        oxygen.tolerant, spores.bool, 
                        Control.Persistence.rank, rCDI.Engraftment.rank, 
                        )) %>% 
  mutate(value_cut = Hmisc::cut2(round(value, 2), cuts = seq(0, 1, 0.2), minmax = T)) %>%
  mutate(facet_title = 'Sp. Features')
  
```

### Bar Plots
Bar plots of per-taxa transmission rates for Control and rCDI data
```{r fig.width=3.5, fig.height=6, message=FALSE, warning=FALSE, paged.print=FALSE}
strain_colors <- c(6, 1, 4, 8, 3, 7, 2)
species_colors <- c(4, 1, 3, 2)
strain_colors_reduced <- c(9, 6, 8, 7, 7 + 5)

p1 <-
  sstr_cases.genus.fraction %>% 
  
  mutate(source = fct_relevel(source, 
                              'Donor', 'Donor Sp.',
                              'Coexistence',
                              'Recipient / Initial Sample','Recipient / Initial Sample Sp.',
                              
                              'Unique to Time Point','Unique to Time Point Sp.',
                              )) %>% 
  mutate(genus = gsub(genus, pattern = '_noname', replacement = '')) %>%  
  
  ggplot(aes(x = fct_reorder(genus, donor_derived.genus, .desc = F), y = f, fill = source)) + 
  geom_bar(stat = 'identity', colour = 'black', size = 0.05) + 
  coord_flip()  +    
  
  guides(fill = guide_legend(ncol = 1)) + 
  facet_grid(cols = vars(Study_Type),
             # rows = vars(control_cluster, rcdi_cluster),
             space = 'free_y',
             scales = 'free') +
  scale_fill_manual(values = colors.discrete[strain_colors]) + 
  tidytext::scale_x_reordered() +
  theme_cowplot(9) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand =c(0,0),
                     labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 9), 
        axis.title = element_blank(),
        axis.text.x = element_text(hjust = 0.5), 
        legend.title=element_blank(),
        legend.position = 'none',
        panel.spacing.x = unit(0.25, "cm")
        # plot.margin = unit(c(0, 0, .9, 0), "cm")
        )
p1
```

Metadata for barplots showing the fraction of species that were found in post-FMT samples within each genus which hold a feature compared to the ones that do not hold a feature. Color coded from white to black (has feature) in 6 steps.
```{r fig.width=7.85, fig.height=6, message=FALSE, warning=FALSE, paged.print=FALSE}
p2 <-
  sstr_cases.genus.fraction.annotated %>% 
  mutate(genus = gsub(genus, pattern = '_noname', replacement = '')) %>% 
  
  mutate(feature = case_when(
  feature == 'gram.bool' ~ 'Gram (+)', 
  feature == 'habit.oral' ~ 'Oral Spp.',
  feature == 'pathogen.human' ~ 'Pathogen',
  feature == 'oxygen.tolerant' ~ 'Oxy. Tolerant',
  feature == 'spores.bool' ~ 'Spore Forming', 
  feature == 'multistrain' ~ 'Multiple Strains', 
  T ~ feature)) %>% 
  
  filter(feature %in% c('Gram (+)','Oral Spp.','Oxy. Tolerant','Spore Forming')) %>% 

  ggplot(aes(y = fct_reorder(genus, donor_derived.genus), 
             x = fct_relevel(feature, 'Gram (+)', 'Spore Forming', 'Oral Spp.', 'Oxy. Tolerant', 
                             ), fill = value_cut)) + 
  geom_tile() +
  coord_equal() + 
  
  facet_grid(cols = vars('')) +
  scale_y_discrete(position = 'right', expand = c(0,0)) +
  
  scale_fill_manual(values = c('white', rev(gray.colors(4)), 'black')) +
  theme_cowplot(9) + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 10), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # 
        axis.title = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.position = 'none',  
        axis.line.y = element_blank())
p2
```


```{r fig.width=7.85, fig.height=6, message=FALSE, warning=FALSE, paged.print=FALSE}
output_name = 'PostSource.Taxa.Genus'
ggsave(p1, device = 'pdf', 
       filename = paste0(results_dir, fig, output_name, '.Data.pdf'), 
       dpi = 300, width = 3.5, height = 6)
ggsave(p2, device = 'pdf', 
       filename = paste0(results_dir, fig, output_name, '.Features.pdf'), 
       dpi = 300, width = 3, height = 6)
```


### Engraftment/Persistence
Engrafted genera in rCDI cases are correlated with persisted taxa in control cases
```{r}
mp_species.rcdi.post.long <- 
  mp_species.rcdi.long %>%
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 
  filter(Sample_Type %in% c('post'), fmt_success) %>% 
  group_by(kingdom, phylum, class, order, family, genus, species, 
    Study, Unique_ID, Days_Since_FMT, Case_Name, Donor.Name, Donor.Unique_ID, Donor.Subject) %>% 
  summarize(rel_abund = sum(rel_abund/100, na.rm = T), .groups = 'drop')
```

Annotate transmission rates with mean post-FMT rel. abund.
```{r fig.height=4, fig.width=4}
genus_transmission_rates.annotated <- 
  mp_species.rcdi.long %>%
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 
  filter(Sample_Type %in% c('post'), fmt_success) %>% 
  group_by(kingdom, phylum, class, order, family, genus,
    Study, Unique_ID, Days_Since_FMT, Case_Name, Donor.Name, Donor.Unique_ID, Donor.Subject) %>% 
  summarize(rel_abund = sum(rel_abund/100, na.rm = T), .groups = 'drop') %>% 
  group_by(kingdom, phylum, class, order, family, genus) %>% 
  summarize(rel_abund.mean = mean(rel_abund, na.rm = T), .groups = 'drop') %>% 
  right_join(., genus_transmission_rates)
```

```{r fig.height=4, fig.width=4}
# correlation metrics
dd <- genus_transmission_rates.annotated %>% 
  filter(!is.na(rel_abund.mean))
dd.cor <- tidy(cor.test(dd$rCDI.Engraftment.rank, dd$Control.Persistence.rank, exact = F, method = 'spearman'))
dd.cor.est <- round(dd.cor$estimate, 2)
dd.cor.sig <- cut(dd.cor$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", ""))

plot <-
  genus_transmission_rates.annotated %>% 
  filter(!is.na(rel_abund.mean)) %>% 
  
  ggplot(aes(rCDI.Engraftment.rank, Control.Persistence.rank, fill = rel_abund.mean)) + 
  geom_point(size = 4, shape = 21) +
  ggrepel::geom_text_repel(aes(label = genus), nudge_x = 0.05,nudge_y = 0.05,
                           data = genus_transmission_rates.annotated %>% 
                    filter((rCDI.Engraftment.rank > 0.01 & n.rCDI.Engraftment > 10 | 
                            Control.Persistence.rank > 0.01 & n.Control.Persistence > 10) & rel_abund.mean > 0.05)
                    ) + 
  coord_fixed() +
  scale_fill_gradient2(name = 'Post-FMT Mean Rel. Abund.',
                       low = colors.discrete[2], mid = colors.discrete[10], high = 'black',
                       midpoint = pseudo_log_trans(1e-6, 10)$transform(0.003),
                       trans = pseudo_log_trans(1e-6, 10), 
                       breaks = c(0.001, 0.01, 0.1),
                       labels = percent_format(), 
                       guide = guide_colorbar(title.position = "top" ,title.hjust = 0.5)) + 
  scale_x_continuous(labels = percent_format_signif) + 
  scale_y_continuous(labels = percent_format_signif) + 
  theme_cowplot() + 
  theme(legend.position = 'top', 
        legend.key.width=unit(1,"cm"),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) + 
  labs(x = 'rCDI - Engraftment rank',
       y = 'Control - Persistence rank') + 
  annotate('text', label = paste0('r=', dd.cor.est,' ', dd.cor.sig), x = .8, y = .05, size = 5)
plot + theme(legend.position = 'none')

legend = cowplot::get_legend(plot)
grid.newpage()
grid.draw(legend)
```

Exporting plot
```{r}
output_name = 'EngraftmentPersistence.Genus'

ggsave(plot + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 4, height = 4,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
ggsave(legend, 
       device = 'pdf', dpi = 300, width = 2.5, height = 1,
       filename = paste0(results_dir, fig, output_name, '.Legend.pdf'))

```

## Per Species

### Format
```{r fig.width=8, fig.height=15, message=FALSE, warning=FALSE, paged.print=FALSE}
sstr_cases.species.fraction <-
  
  sstr_cases.species %>% 
  
  ## filter genera
  filter(species %in% species_importance$species) %>% 
  pivot_longer(names_to = 'Study_Type', values_to = 'n', cols = c('Control','rCDI')) %>% 

  ## calculate fraction of events
  group_by(Study_Type, species, source.simple) %>%
  mutate(source.n = sum(n, na.rm = T)) %>% 
  
  group_by(Study_Type, species) %>%
  mutate(total = sum(n, na.rm = T), 
         f = n / total,
         source.f = source.n / total) %>% 
  
  ## get species-level donor-derived for ordering
  group_by(species) %>% 
  mutate(donor_derived.species = max(ifelse(Study_Type == 'rCDI' & 
                                          source.simple == 'Donor', 
                                          source.f, 0), na.rm = T)) %>% 
  ungroup()  


sstr_cases.species.fraction.annotated <- 
  sstr_cases.species.fraction %>% 
  
  # add species metrics for present species
  left_join(., species_feature_rates) %>%
  left_join(., species_transmission_rates) %>%
  
  pivot_longer(names_to = 'feature', values_to = 'value', 
               cols = c(gram.bool, habit.oral, 
                        oxygen.tolerant, spores.bool, 
                        Control.Persistence.rank, rCDI.Engraftment.rank, 
                        )) %>% 
  mutate(value_cut = Hmisc::cut2(round(value, 2), cuts = seq(0, 1, 0.2), minmax = T)) %>%
  mutate(facet_title = 'Sp. Features')
```

### Bar Plots
Bar plots of per-taxa transmission rates for Control and rCDI data
```{r fig.width=5, fig.height=15, message=FALSE, warning=FALSE, paged.print=FALSE}
strain_colors <- c(6, 1, 4, 8, 3, 7, 2)
species_colors <- c(4, 1, 3, 2)
strain_colors_reduced <- c(9, 6, 8, 7, 7 + 5)

p1 <-
  sstr_cases.species.fraction %>% 
  
  mutate(source = fct_relevel(source, 
                              'Donor', 'Donor Sp.',
                              'Coexistence',
                              'Recipient / Initial Sample','Recipient / Initial Sample Sp.',
                              
                              'Unique to Time Point','Unique to Time Point Sp.',
                              )) %>% 
  mutate(species = gsub(species, pattern = '_noname', replacement = '')) %>%  
  
  ggplot(aes(x = fct_reorder(species, donor_derived.species, .desc = F), y = f, fill = source)) + 
  geom_bar(stat = 'identity', colour = 'black', size = 0.05) + 
  coord_flip()  +    
  
  guides(fill = guide_legend(ncol = 1)) + 
  facet_grid(cols = vars(Study_Type),
             # rows = vars(control_cluster, rcdi_cluster),
             space = 'free_y',
             scales = 'free') +
  scale_fill_manual(values = colors.discrete[strain_colors]) + 
  tidytext::scale_x_reordered() +
  theme_cowplot(9) +
  scale_x_discrete(expand = c(0,0)) +
  scale_y_continuous(expand =c(0,0),
                     labels = scales::percent_format(), breaks = c(0, 0.5, 1)) +
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 9), 
        axis.title = element_blank(),
        axis.text.x = element_text(hjust = 0.5), 
        legend.title=element_blank(),
        legend.position = 'none',
        panel.spacing.x = unit(0.25, "cm")
        # plot.margin = unit(c(0, 0, .9, 0), "cm")
        )
p1
```

Metadata for barplots showing the fraction of species that were found in post-FMT samples within each species which hold a feature compared to the ones that do not hold a feature. Color coded from white to black (has feature) in 6 steps.
```{r fig.width=3, fig.height=15, message=FALSE, warning=FALSE, paged.print=FALSE}
p2 <-
  sstr_cases.species.fraction.annotated %>% 
  mutate(species = gsub(species, pattern = '_noname', replacement = '')) %>% 
  
  mutate(feature = case_when(
  feature == 'gram.bool' ~ 'Gram (+)', 
  feature == 'habit.oral' ~ 'Oral Spp.',
  feature == 'pathogen.human' ~ 'Pathogen',
  feature == 'oxygen.tolerant' ~ 'Oxy. Tolerant',
  feature == 'spores.bool' ~ 'Spore Forming', 
  feature == 'multistrain' ~ 'Multiple Strains', 
  T ~ feature)) %>% 
  
  filter(feature %in% c('Gram (+)','Oral Spp.','Oxy. Tolerant','Spore Forming')) %>% 

  ggplot(aes(y = fct_reorder(species, donor_derived.species), 
             x = fct_relevel(feature, 'Gram (+)', 'Spore Forming', 'Oral Spp.', 'Oxy. Tolerant', 
                             ), fill = value_cut)) + 
  geom_tile() +
  coord_equal() + 
  
  facet_grid(cols = vars('')) +
  scale_y_discrete(position = 'right', expand = c(0,0)) +
  
  scale_fill_manual(values = c('white', 'black')) +
  theme_cowplot(9) + 
  theme(strip.background = element_blank(), 
        strip.text = element_text(size = 10), 
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1), # 
        axis.title = element_blank(),
        axis.ticks.y = element_blank(), 
        legend.position = 'none',  
        axis.line.y = element_blank())
p2
```


```{r fig.width=7.85, fig.height=6, message=FALSE, warning=FALSE, paged.print=FALSE}
output_name = 'PostSource.Taxa.Species'
ggsave(p1, device = 'pdf', 
       filename = paste0(results_dir, fig, output_name, '.Data.pdf'), 
       dpi = 300, width = 5, height = 15)
ggsave(p2, device = 'pdf', 
       filename = paste0(results_dir, fig, output_name, '.Features.pdf'), 
       dpi = 300, width = 3, height = 15)
```


### Engraftment/Persistence
Engrafted species in rCDI cases are correlated with persisted taxa in control cases
```{r}
mp_species.rcdi.post.long <- 
  mp_species.rcdi.long %>%
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 
  filter(Sample_Type %in% c('post'), fmt_success) %>% 
  group_by(kingdom, phylum, class, order, family, genus, species, 
    Study, Unique_ID, Days_Since_FMT, Case_Name, Donor.Name, Donor.Unique_ID, Donor.Subject) %>% 
  summarize(rel_abund = sum(rel_abund/100, na.rm = T), .groups = 'drop')
```

Annotate transmission rates with mean post-FMT rel. abund.
```{r fig.height=4, fig.width=4}
species_transmission_rates.annotated <- 
  mp_species.rcdi.long %>%
  filter(kingdom == 'Bacteria', !grepl(species, pattern = 'unclassified')) %>% 
  filter(Sample_Type %in% c('post'), fmt_success) %>% 
  group_by(kingdom, phylum, class, order, family, genus, species) %>% 
  summarize(rel_abund.mean = mean(rel_abund/100, na.rm = T), .groups = 'drop') %>% 
  right_join(., species_transmission_rates)
```

```{r fig.height=4, fig.width=4}
# correlation metrics
dd <- species_transmission_rates.annotated %>% 
  filter(!is.na(rel_abund.mean))
dd.cor <- tidy(cor.test(dd$rCDI.Engraftment.rank, dd$Control.Persistence.rank, exact = F, method = 'spearman'))
dd.cor.est <- round(dd.cor$estimate, 2)
dd.cor.sig <- cut(dd.cor$p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                    labels = c("***", "**", "*", ""))

plot <-
  species_transmission_rates.annotated %>% 
  filter(!is.na(rel_abund.mean)) %>% 
  
  ggplot(aes(rCDI.Engraftment.rank, Control.Persistence.rank, fill = rel_abund.mean)) + 
  geom_point(size = 4, shape = 21) +
  ggrepel::geom_text_repel(aes(label = species), data = species_transmission_rates.annotated %>% 
                    filter((rCDI.Engraftment.rank > 0.01 & n.rCDI.Engraftment >= 10 | 
                            Control.Persistence.rank > 0.01 & n.Control.Persistence >= 10) & rel_abund.mean > 0.03)
                    ) + 
  coord_fixed() +
  scale_fill_gradient2(name = 'Post-FMT Mean Rel. Abund.',
                       low = colors.discrete[2], mid = colors.discrete[10], high = 'black',
                       midpoint = pseudo_log_trans(1e-6, 10)$transform(0.001),
                       trans = pseudo_log_trans(1e-6, 10), 
                       breaks = c(0.0001, 0.001, 0.01),
                       labels = percent_format(), 
                       guide = guide_colorbar(title.position = "top" ,title.hjust = 0.5)) + 
  scale_x_continuous(labels = percent_format_signif) + 
  scale_y_continuous(labels = percent_format_signif) + 
  theme_cowplot() + 
  theme(legend.position = 'top', 
        legend.key.width=unit(1,"cm"),
        legend.key.size = unit(0.35, "cm"),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10)) + 
  labs(x = 'rCDI - Engraftment rank',
       y = 'Control - Persistence rank') + 
  annotate('text', label = paste0('r=', dd.cor.est,' ', dd.cor.sig), x = 0, y = 1, hjust = 0, size = 5)
plot + theme(legend.position = 'none')

legend = cowplot::get_legend(plot)
grid.newpage()
grid.draw(legend)
```

Exporting plot
```{r}
output_name = 'EngraftmentPersistence.Species'

ggsave(plot + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 4, height = 4,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
ggsave(legend, 
       device = 'pdf', dpi = 300, width = 2.5, height = 1,
       filename = paste0(results_dir, fig, output_name, '.Legend.pdf'))

```

# Post-FMT Species Rel. Abund.
Rel. abund. of post-FMT species is dependent on which strain survived the FMT treatment.

## Format
```{r}
competing_cases <-
  sstr_cases %>% 
  filter(Study_Type == 'rCDI') %>% 
  filter(replace_na(rel_abund.recipient, 0) > 0 & replace_na(rel_abund.donor, 0) > 0) %>% 
  select(taxonomy.levels, event, source, analysis_level, combined_resolution,
         starts_with('existed'),
         Study, Case_Name, Days_Since_FMT.post, Study_Type, fmt_success, 
         rel_abund.recipient, rel_abund.post, rel_abund.donor, 
         ends_with('.mvs'), ends_with('.overlap'),
         starts_with('f_poly.'), starts_with('n_covered.')) %>% 
         left_join(., microbe.metadata)
```

Stats
```{r}
competing_cases %>% 
  filter(analysis_level != 'species') %>% 
  group_by(source) %>% 
  tally() %>% 
  ungroup() %>% 
  mutate(f = n / sum(n, na.rm = T)) %>% 
  column_to_rownames('source') %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column('metric')
```
Replacement = 207 (50.7%)
Resistence = 119 (29.2%)
New/Unique = 57 (13.9%)
Co-Existence = 25 (6.1%)
Total = 207 + 119 + 57 + 25 = 408

Export Table
```{r}
output_name = 'Competing_Cases'
competing_cases %>%
  filter(analysis_level != 'species') %>% 
  select(Study, Case_Name, Days_Since_FMT.post, fmt_success, 
         species, source, starts_with('rel_abund'), 
         ends_with('mvs'), ends_with('overlap'), starts_with('n_covered')) %>% 
  write_tsv(., paste0(results_dir, fig, output_name, '.tsv'))
```

## Strain-only
Events with all three R/D/P comparisons at strain-level
```{r fig.width=4,fig.height=8}
bb = c(1e-5, 1e-3, 1e-1)

competing_cases.strain <-
  competing_cases %>%
  filter(replace_na(rel_abund.post, 0) > 0) %>% 
  filter(n_covered.donor > min_overlap & 
         n_covered.recipient > min_overlap & 
         n_covered.post > min_overlap) %>% 
  filter(source %in% c('donor', 'self', 'both'), source != 'same_species') %>%
  pivot_longer(names_to = 'rel_abund.source', values_to = 'rel_abund', 
               cols = c('rel_abund.recipient', 'rel_abund.donor')) %>% 
  mutate(rel_abund.source = case_when(
    rel_abund.source == 'rel_abund.recipient' ~ 'Recipient', 
    rel_abund.source == 'rel_abund.donor' ~ 'Donor',
    T ~ rel_abund.source)) %>% 
  mutate(rel_abund.source = fct_relevel(as.factor(rel_abund.source), 'Recipient')) %>% 
  mutate(source = case_when(
    source == 'both' ~ 'Coexistence',
    source == 'self' ~ 'Recipient-Specific', 
    source == 'donor' ~ 'Donor-Specific',
    T ~ source
  ),
  source = fct_relevel(source, 'Recipient-Specific', 'Donor-Specific')) 

min_cutoff = min(c(min(competing_cases.strain$rel_abund, na.rm = T), 
                   min(competing_cases.strain$rel_abund.post, na.rm = T))) / 100
```


The relative abundance of the post-FMT species correlates with the relative abundance of the recipient but not donor species, if the recipient species persisted and vice versa.
```{r fig.width=4,fig.height=8}
plot <- 
  competing_cases.strain %>% 
  ggplot(aes((rel_abund.post / 100), (rel_abund / 100))) +
  geom_point(aes(fill = rel_abund.source), shape = 21, size = 3, col = 'black', alpha = 1) +
  scale_x_log10(label = percent_format_signif,
                breaks = bb,
                expand = c(0.1,0),
                limits = c(min_cutoff, 1)
                ) +
  scale_y_log10(label = percent_format_signif,
                breaks = bb,
                expand = c(0.1,0),
                limits = c(min_cutoff, 200) # space for stat_cor
                ) +
  theme_cowplot() +
  stat_cor(method = 'spearman', cor.coef.name = 'r',
           aes(col = rel_abund.source, 
               label = paste0(str_replace(str_remove_all(str_to_lower(..r.label..), 
                                                         pattern = ' '), 
                                          pattern = ',', replacement = ', '))),
               output.type="text", show.legend = F) +
  geom_smooth(method = 'lm', se =F, aes(col = rel_abund.source), show.legend = F) +
  scale_color_manual(values = colors.discrete[c(3, 1)]) +
  scale_fill_manual(values = colors.discrete[c(3, 1)],
                    guide = guide_legend(direction = "horizontal", 
                             title.position = "top", 
                             title = 'Rel. Abund.', 
                             title.theme = element_text(angle = 270, size = 14),
                             
                             label.position="bottom", 
                             label.hjust = 0, 
                             label.vjust = 0.5,
                             label.theme = element_text(angle = 270)
                             )
                    ) +
  labs(x = 'Post-FMT') + 
  theme(aspect.ratio = 1,
        legend.position = 'left',
        strip.background = element_blank(), 
        axis.text.x = element_text(hjust = .65), 
        axis.title.y = element_blank()) +
  facet_wrap(~source, ncol = 1) +
  geom_abline(intercept = 0, slope = 1, linetype= 'dashed')

plot + theme(legend.position = 'none')

legend = cowplot::get_legend(plot)
grid.newpage()
grid.draw(legend)
```

## Combined, 2/3
Events with 2/3 R/D/P comparisons at strain-level
```{r fig.width=8,fig.height=4}
bb = c(1e-5, 1e-3, 1e-1)

competing_cases.combined <-
  competing_cases %>%
  filter(replace_na(rel_abund.post, 0) > 0) %>% 
  filter(source %in% c('donor', 'self', 'both'), source != 'same_species') %>%
  pivot_longer(names_to = 'rel_abund.source', 
               values_to = 'rel_abund', 
               cols = c('rel_abund.recipient', 'rel_abund.donor')) %>% 
  mutate(rel_abund.source = case_when(
    rel_abund.source == 'rel_abund.recipient' ~ 'Recipient', 
    rel_abund.source == 'rel_abund.donor' ~ 'Donor',
    T ~ rel_abund.source)) %>% 
  mutate(rel_abund.source = fct_relevel(as.factor(rel_abund.source), 'Recipient')) %>% 
  mutate(source = case_when(
    source == 'both' ~ 'Coexistence',
    source == 'self' ~ 'Recipient-Derived', 
    source == 'donor' ~ 'Donor-Derived',
    T ~ source
  ),
  source = fct_relevel(source, 'Recipient-Derived', 'Donor-Derived')) 

min_cutoff = min(c(min(competing_cases.combined$rel_abund, na.rm = T), 
                   min(competing_cases.combined$rel_abund.post, na.rm = T))) / 100
```


```{r fig.width=4,fig.height=8}
plot <-
  competing_cases.combined %>% 
  ggplot(aes((rel_abund.post / 100), (rel_abund / 100))) +
  geom_point(aes(fill = rel_abund.source), shape = 21, size = 3, col = 'black', alpha = 1) +
  scale_x_log10(label = percent_format_signif,
                breaks = bb,
                expand = c(0.1,0),
                limits = c(min_cutoff, 1)
                ) +
  scale_y_log10(label = percent_format_signif,
                breaks = bb,
                expand = c(0.1,0),
                limits = c(min_cutoff, 200) # space for stat_cor
                ) +
  theme_cowplot() +
  stat_cor(method = 'spearman', 
           aes(col = rel_abund.source, 
               label = paste0(str_replace(str_remove_all(str_to_lower(..r.label..), 
                                                         pattern = ' '), 
                                          pattern = ',', replacement = ', '))),
               output.type="text", show.legend = F) +
  geom_smooth(method = 'lm', se =F, aes(col = rel_abund.source), show.legend = F) +
  scale_color_manual(values = colors.discrete[c(3, 1)]) +
  scale_fill_manual(values = colors.discrete[c(3, 1)],
                    guide = guide_legend(direction = "horizontal", 
                             title.position = "top", 
                             title = 'Rel. Abund.', 
                             title.theme = element_text(angle = 270, size = 14),
                             
                             label.position="bottom", 
                             label.hjust = 0, 
                             label.vjust = 0.5,
                             label.theme = element_text(angle = 270)
                             )
                    ) +
  labs(x = 'Post-FMT') + 
  theme(aspect.ratio = 1,
        legend.position = 'left',
        strip.background = element_blank(), 
        axis.text.x = element_text(hjust = .65), 
        axis.title.y = element_blank()) +
  facet_wrap(~source, ncol = 1) +
  geom_abline(intercept = 0, slope = 1, linetype= 'dashed')

plot + theme(legend.position = 'none')

legend = cowplot::get_legend(plot)
grid.newpage()
grid.draw(legend)
```

Exporting plot
```{r}
output_name = 'StrainSource.RelAbund.rCDI'


ggsave(plot + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 3.5, height = 7,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
ggsave(legend, 
       device = 'pdf', dpi = 300, width = 1, height = 3,
       filename = paste0(results_dir, fig, output_name, '.Legend.pdf'))
```

# Post-FMT Presence/Absence (Strain)
Using a random forest classifier to predict whether strains found in the donor will be found post-FMT in the patient

## Format
```{r}
scale_gelman <- function(x) {
  return((x - mean(x, na.rm = T)) / (2 * sd(x, na.rm = T)))
}

sstr_cases.format <- 
  sstr_cases %>%

  filter(Study_Type %in% c('rCDI') ) %>% 
  filter(kingdom == 'Bacteria') %>% 

  mutate(source = ifelse(grepl(species, pattern = 'unclassified'),  'Unclassified', source)) %>%
  mutate(source = case_when(
    analysis_level == 'species' & source == 'self' ~ 'Self Sp.',
    analysis_level == 'species' & source == 'donor' ~ 'Donor Sp.',
    analysis_level == 'species' & source == 'unique' ~ 'Unique Sp.',
    T ~ source
  )) %>% 
  
  left_join(., microbe.metadata) %>% 
  
  mutate_at(.vars = vars(rel_abund.donor, rel_abund.post, rel_abund.recipient),
            .funs = funs(. / 100)) %>% 
  mutate(Donor = replace_na(rel_abund.donor > 0, F), 
         PostTreatment = replace_na(rel_abund.post > 0, F), 
         Pre = replace_na(rel_abund.recipient > 0, F),
         Both = Donor & Pre, 
         Any = Donor | Pre) %>% 
  
  # filter only species that existed in donor or both
  filter(!grepl(species, pattern = 'unclassified'), species %in% microbe.metadata$species) %>% 
  filter(Donor | Both) %>%  # n_covered.donor > min_overlap

  mutate(Engrafted = case_when(
    
    # shared_strain
    `Donor/Post-FMT.mvs` > min_similarity &
      `Donor/Post-FMT.overlap` > min_overlap ~ T,
    
    # other_strain
    `Donor/Post-FMT.mvs` < min_similarity &
      `Donor/Post-FMT.overlap` > min_overlap ~ F,
    
    # unique_species
    replace_na(`Donor/Post-FMT.overlap`, 0) < min_overlap &
      !PostTreatment ~ F
    )) %>% 
  
  # filter only events for which we know outcome
  filter(!is.na(Engrafted)) %>% 
  
    
  rename(OxyTol = 'oxygen.tolerant', 
         OxyClass = 'oxygen.class',
         DaysSinceFMT = 'Days_Since_FMT.post', 
         Case = 'Case_Name',
         AbundanceDonor = 'rel_abund.donor',
         AbundanceRecipient = 'rel_abund.recipient',
         AbundancePost = 'rel_abund.post') %>% 
  
  mutate(Habitat = ifelse(habit.oral, 'Oral', 'Not-Oral'),
         OxyTol = ifelse(OxyTol, 'Tolerant', 'Not-Tolerant')) %>% 
  mutate(Detected = case_when(
    Pre & !Donor ~ 'Recipient', 
    !Pre & Donor ~ 'Donor',
    Pre & Donor ~ 'Both')) %>% 
  mutate(Specificity = ifelse(Detected == 'Both', F, T)) %>% 
  mutate_at(.vars = vars(Case, 
                         Specificity, Study, 
                         Detected, Engrafted,
                         phylum, class, order, family, 
                         OxyTol, OxyClass, Habitat, PostTreatment,
                         habit.site, gram_stain, spore_forming
                         ),
            .funs = funs(as.factor(.))) %>% 
  mutate(
  Detected = fct_relevel(Detected, 'Recipient'),
  OxyClass = fct_relevel(OxyClass, 'aerobe'),
  Habitat = fct_relevel(Habitat, 'Not-Oral'),
  OxyTol = fct_relevel(OxyTol, 'Not-Tolerant'),
  Phylum = fct_relevel(phylum, 'Firmicutes')
  ) %>% 
  
  mutate(OxyTol = as.factor(replace_na(as.character(OxyTol), 'Unknown')),
         OxyClass = as.factor(replace_na(as.character(OxyClass), 'unknown'))) %>% 
  
  select(phylum, class, order, family, genus, species, Study, Engrafted, PostTreatment, AbundancePost, Detected, AbundanceDonor, AbundanceRecipient, Habitat, OxyTol, DaysSinceFMT, Case, OxyClass, habit.site, gram_stain, spore_forming)
```

Focussing on donor-derived taxa
```{r}
d <- 
  sstr_cases.format %>% 
  
  mutate_at(.vars = vars(AbundanceDonor, AbundanceRecipient),
            .funs = funs(pseudo_log_trans(1e-6, 10)$transform(.))) %>% 
  mutate(AbundanceRatio = AbundanceDonor - AbundanceRecipient, # diff is ratio at log
         DaysSinceFMT = log10(DaysSinceFMT)
         )  %>% 
   
  mutate_at(.vars = vars(starts_with('Abundance'), DaysSinceFMT),
            .funs = funs(scale_gelman(.))) %>% 

  select(-AbundancePost, 
         -PostTreatment, -Case, 
         -phylum, 
         -class,
         -genus,
         -species,
         -habit.site
         )
```

## Random Forest Model
Randomly sample 80/20 split, use RF to model engraftment, test on 20% hold-out data
```{r}
seed = 30
set.seed(seed) 
trainIndex <- createDataPartition(d$Engrafted, p = .8, 
                                  list = FALSE, 
                                  times = 1)[,1]
train <- d[ trainIndex,]
test <- d[-trainIndex,]
```

```{r}
unlist(lapply(list(summary(train$Engrafted)), function(x) x / sum(x)))
nrow(train)
unlist(lapply(list(summary(test$Engrafted)), function(x) x / sum(x)))
nrow(test)
```


Model
```{r}
set.seed(seed)

rf <- randomForest(
  Engrafted ~ .,
  na.action = na.roughfix,
  data = train, 
  ntree = 500, nodesize = 1,
  mtry = 3, importance = T
)

# imp based on validation data
imp <-
  randomForest::importance(rf, type = 1, scale = F, )
```

## Performance
```{r fig.width=3.5, fig.height=3.5}
options(yardstick.event_first = FALSE)

test$rf <- predict(rf, test, type = 'prob')[,2]

roc <-
  test %>% 
  ggplot(aes(m = rf, d = Engrafted)) + 
  geom_roc(n.cuts = 0, labels = F)

auPR <- 
  test %>% 
  pr_auc(data = ., truth = Engrafted, rf, na_rm = T) %>% 
  pull(.estimate) %>% round(., digits = 3)

test %>% 
  pr_curve(data = ., truth = as.factor(Engrafted), rf) %>% 
  ggplot(aes(y = precision, x = recall)) + 
  geom_line(size = 1) + 
  coord_fixed() + 
  theme_cowplot() +
  scale_color_manual(values = colors.discrete) + 
  scale_x_continuous(labels = percent_format(), limits = c(0, 1)) + 
  scale_y_continuous(labels = percent_format(), limits = c(0, 1), sec.axis = dup_axis()) + 
  labs(x = 'Recall', y = 'Precision') + 
  geom_abline(intercept = 1, slope = -1, linetype = 'dashed') + 
  theme(legend.position = 'none',
        axis.line.y.right = element_blank(),
        axis.ticks.y.right = element_blank(),
        axis.text.y.right = element_blank(),
        axis.title.y.left = element_blank())

auROC <- calc_auc(roc)$AUC %>% round(., digits = 3)

auROC
auPR

plot <- 
  roc + 
  annotate('text', label = paste0('auROC: ', auROC), x = 0.4, y = 0.15, hjust = 0, size = 5) + 
  annotate('text', label = paste0('auPR: ', auPR), x = 0.4, y = 0.05, hjust = 0, size = 5) + 
  theme_cowplot() +
  coord_fixed() + 
  scale_x_continuous(labels = percent_format(), limits = c(0, 1)) +
  scale_y_continuous(labels = percent_format(), limits = c(0, 1)) +
  geom_abline(slope = 1, linetype = 'dashed') + 
  labs(x = 'False Positive', y = 'True Positive') + 
  theme(aspect.ratio = 1)
plot 
```

Exporting plot
```{r}
output_name = 'PostPresence.Strain.RF.ROC'

ggsave(plot + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 4, height = 4,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
```


## Variable Importance
Permuted vip
```{r}
rf.vip_acc <- as.data.frame(rf$importance) %>% 
  rownames_to_column('feature') %>% 
  select(feature, MeanDecreaseAccuracy)
rf.vip_sd <- as.data.frame(rf$importanceSD) %>% 
  rownames_to_column('feature') %>% 
  select(feature, MeanDecreaseAccuracy) %>% 
  mutate(sd = MeanDecreaseAccuracy) %>% 
  select(feature, sd)

rf.vip <- full_join(rf.vip_acc, rf.vip_sd, by = 'feature')
```


```{r}
feat_order <- c('Taxonomy', 'Species\nAbundance', 'Clinical', 'Physiology')
plot <- rf.vip %>% 
  rename(importance = 'MeanDecreaseAccuracy') %>% 
  mutate(feature = as.character(feature)) %>% 
  mutate(category = case_when(
    feature %in% c('kingdom', 'phylum', 'class', 'order', 'family') ~ 'Taxonomy', 
    feature %in% c('DaysSinceFMT', 'Study') ~ 'Clinical', 
    feature %in% c('AbundanceDonor', 'AbundanceRecipient', 'AbundanceRatio', 'Detected') ~ 'Species\nAbundance', 
    feature %in% c('Habitat', 'OxyClass', 'OxyTol', 'gram_stain', 'spore_forming') ~ 'Physiology', 
    T ~ 'other'
  )) %>% 
  
  mutate(feature = case_when(
    grepl(feature, pattern = '^Abundance') ~ str_split_fixed(feature, 'Abundance', 2)[,2], 
    feature == 'DaysSinceFMT' ~ 'Days Since FMT', 
    feature == 'Detected' ~ 'Donor-Specificity', 
    feature == 'GCcontent' ~ '%GC-Content', 
    feature == 'OxyTol' ~ 'Oxygen Tolerance', 
    feature == 'OxyClass' ~ 'Oxygen Usage', 
    feature == 'Habitat' ~ 'Oral Habitat',
    feature == 'family' ~ 'Family',
    feature == 'order' ~ 'Order',
    feature == 'gram_stain' ~ 'Gram Staining',
    feature == 'spore_forming' ~ 'Spore Formation',
    T ~ feature
  )) %>% 
  
  ggplot(aes(fct_reorder(feature, importance), importance)) + 
  geom_bar(stat = 'identity', fill = 'black', width = 0.1) +
  geom_point(size = 3) + 
  coord_flip() + 
  scale_y_continuous(labels = percent_format_signif) + 
  facet_grid(rows = vars(fct_relevel(category, feat_order)), scales = 'free_y', space = 'free_y') +
  labs(x= '', y = 'Mean Decrease in Accuracy') +
  theme_cowplot() + 
  theme(strip.background = element_blank(), 
        strip.text.y = element_text(angle = 0),
        # strip.text.y = element_blank()
        )
  plot
```

Exporting plot
```{r}
output_name = 'PostPresence.Strain.RF.VarImp'

ggsave(plot + theme(legend.position = 'none'), 
       device = 'pdf', dpi = 300, width = 5, height = 4,
       filename = paste0(results_dir, fig, output_name, '.Plot.pdf'))
```


## Abundance in relation to Engraftment Frequency

Format
```{r}
sstr_cases.select <- 
  sstr_cases.format %>% 

  select(species, Case, AbundanceDonor, AbundanceRecipient, Engrafted) %>% 
  
  mutate(AbundanceRatio = 
           pseudo_log_trans(1e-6, 10)$transform(AbundanceDonor) - 
           pseudo_log_trans(1e-6, 10)$transform(AbundanceRecipient)) %>% 
  select(-species, -Case)
```

What fraction of events are engrafted, given that the relative abundance of the donor strain is >=x%?
```{r fig.width=3.5, fig.height=3.5}
sstr_cases.select.AbundanceDonor <- 
  sstr_cases.select %>% arrange(desc(AbundanceDonor)) %>% 
  select(AbundanceDonor, Engrafted)

get_donor_fraction_by_abundance <- function(min_abund) {
  x <- sstr_cases.select.AbundanceDonor %>% 
    filter(AbundanceDonor >= min_abund)
  return(sum(x$Engrafted == T) / nrow(x))
}

r.abund <- c()
abundance_steps <- 10^(seq(
  log10(min(sstr_cases.select.AbundanceDonor$AbundanceDonor)), 
  log10(sstr_cases.select.AbundanceDonor$AbundanceDonor[15]), # at least 15 observations
  length.out = 100))
for (as in abundance_steps) {r.abund <- c(r.abund, get_donor_fraction_by_abundance(as))}

r.abund.df <- 
data.frame(min_rel_abund = abundance_steps, 
           engrafted_fraction = r.abund)
plot.donor <- 
  r.abund.df %>% 
  ggplot(aes(x = min_rel_abund, 
             y = engrafted_fraction)) + 
  geom_line() +
  scale_x_continuous(trans = 'log10', 
                     labels = c('', '.001%', '', '.1%', '', '10%', ''),
                     breaks = c(1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1)
                     ) + 
  scale_y_continuous(labels = percent_format_signif) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1) + 
  labs(x = 'Rel. Abund. in Donor', y = 'Cumulative\nEngraftment Frequency')
plot.donor
```

What fraction of events are engrafted, given that the ratio of relative abundances of the donor to recipient is >=x?
```{r fig.width=3.5, fig.height=3.5}
sstr_cases.select.AbundanceRatio <- 
  sstr_cases.select %>% arrange(desc(AbundanceRatio)) %>% 
  filter(AbundanceRecipient > 0) %>% 
  select(AbundanceRatio, Engrafted)

get_donor_fraction_by_ratio <- function(min_ratio) {
  x <- sstr_cases.select.AbundanceRatio %>% filter(AbundanceRatio >= min_ratio)
  return(sum(x$Engrafted == T) / nrow(x))
}


r.ratio <- c()
abundance_steps <- seq(
  min(sstr_cases.select.AbundanceRatio$AbundanceRatio), 
  sstr_cases.select.AbundanceRatio$AbundanceRatio[15], # at least 15 observations
  length.out = 100)
for (as in abundance_steps) {r.ratio <- c(r.ratio, get_donor_fraction_by_ratio(as))}

r.ratio.df <- 
data.frame(min_ratio = abundance_steps, 
           engrafted_fraction = r.ratio)

plot.ratio <-
  r.ratio.df %>% 
  ggplot(aes(x = min_ratio, 
             y = engrafted_fraction)) + 
  geom_line() +
  scale_y_continuous(labels = percent_format_signif) + 
  theme_cowplot() + 
  theme(aspect.ratio = 1) + 
  labs(x = 'Donor/Recipient (Log-Ratio)', y = 'Cumulative\nEngraftment Frequency')
plot.ratio
```


Exporting plots
```{r}
output_name = 'EngraftmentFrequency'

ggsave(plot.donor, 
       device = 'pdf', dpi = 300, width = 4, height = 4,
       filename = paste0(results_dir, fig, output_name, '.Donor.Plot.pdf'))
ggsave(plot.ratio, 
       device = 'pdf', dpi = 300, width = 4, height = 4,
       filename = paste0(results_dir, fig, output_name, '.Ratio.Plot.pdf'))
```
