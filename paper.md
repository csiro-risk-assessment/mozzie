---
title: "Mozzie: a computationally efficient simulator for the spatio-temporal modelling
  of mosquitoes"
tags:
- Python
- Cython
- malaria
- mosquito
- risk analysis
date: "1 July 2024"
output: pdf_document
authors:
- name: Andy Wilkins
  orcid: "0000-0001-6472-9560"
  equal-contrib: yes
  affiliation: '1'
- name: Nicholas J. Beeton
  orcid: "0000-0002-8513-3165"
  equal-contrib: yes
  affiliation: '1'
- name: "Maud El-Hachem"
  orcid: "0000-0002-9312-6109"
  equal-contrib: no
  affiliation: '1'
- name: Keith R Hayes
  orcid: "0000-0003-4094-3575"
  equal-contrib: no
  affiliation: '1'
- name: Geoffrey R. Hosack
  orcid: "0000-0002-6462-6817"
  equal-contrib: no
  affiliation: '1'
bibliography: paper.bib
affiliations:
- name: Commonwealth Scientific and Industrial Organisation (CSIRO), Australia
  index: 1
---

<!---
Process locally by cloning whedon (git clone git@github.com:openjournals/whedon.git) and then:
pandoc --citeproc --csl ../whedon/resources/apa.csl --pdf-engine=xelatex --variable colorlinks=true paper.md -o paper.pdf
-->

<!---
See [instructions](https://joss.readthedocs.io/en/latest/paper.html) for instructions about formulae, references, tables, figures, markdown, footnotes, headings, etc.  See [example](https://joss.readthedocs.io/en/latest/example_paper.html) for an example.
-->


# Summary

`Mozzie` enables simulation of the lifecycle and spatial spread of mosquitoes.  `Mozzie` can be used to assess risks associated with disease-control strategies at local, regional or continental scales.  Most particularly, strategies involving genetic alterations of mosquitoes to eliminate malaria, are of prime interest.

More technically, `Mozzie` simulates a population-dynamics model that uses differential equations or delay differential equations [@BOHNER2018114;@El-Hachem2024] to describe the spread and persistence of mosquitoes that may be genetically altered.  Genetic alterations are flexibly modelled: these can involve any number of alleles; Mendelian or non-Mendelian inheritance, including gene drives; they can be self-limiting or self-sustaining; and can include the emergence of resistant allelles.  The model allows simulation of $N$ mosquito species.   It incorporates mate-choice, hybridisation and intra-specific competition that occur within complexes of mosquito species [@BEETON2020110072].  This fills a gap that currently exists among similar models, allowing researchers to assess potential transfer of the genetic alterations between (sub-)species.

`Mozzie` supports spatial and temporal variations in lifecyle parameters, and local diffusion and wind-assisted, long range, advection.  For example, wind patterns and the capacity of the landscape to support mosquitoes can vary spatially and temporally, reflecting daily variations, seasonality, and local conditions.

Conversely, `Mozzie` does not contain human agents, nor does it consider the effect of genetic control strategies on the prevalence of pathogens such as the malaria parasite, among human or animal populations.

`Mozzie` has been used by the authors to simulate the spread across sub-Saharan Africa of a theoretical, population-modifying, gene drive in _Anopheles gambiae s.s._ and _Anopheles coluzzii_ [@beetonplos] (that paper also describes the mathematics of a particular mosquito lifecycle model that is contained in `Mozzie`).  It has also been used to predict the spread of Target Malaria's Paternal Male Bias construct [@Galizi2014] following a proposed field-release of genetically modified _Anopheles coluzzii_ male mosquitoes in Burkina Faso [@Hosack2023].


<!---
AN ALTERATE PARAGRAPH:

`Mozzie` enables simulation of the lifecycle and spatial spread of mosquitoes.  Mosquitoes transmit malaria, and `Mozzie` is designed to support risk assessments associated with genetic alterations of mosquitoes aimed at eliminating malaria.  Treatment of the mosquito lifecycle includes competition and breeding between multiple (sub-)species, allowing researchers to assess potential transfer of the genetic alteration to other (sub-)species.  To assess the spatial spread of wild-type and genetically-altered mosquitoes, `Mozzie` can simulate the local diffusion of individuals through the landscape, and long-distance dispersal via wind.  `Mozzie` can be used for risk assessments at the local, regional or continental scale.  For example, wind patterns and the capacity of the landscape to support mosquitoes can vary spatially and temporally, reflecting seasonality, local conditions and/or daily variations, etc.  Genetic alterations can be flexibly modelled with any number of alleles and the capacity for Mendelian or non-Mendelian inheritance, including gene drives.
-->

# State of the field

Alternatives to `Mozzie` include:

- MIT's [HYDREMATS](https://eltahir.mit.edu/models/hydremats/) software [@bomblies2009].  HYDREMATS is a coupled hydrology and entomology model that uses an agent-based approach for mosquito-human dynamics, and focuses on high-resolution village-scale understanding of malaria without genetically-modified mosquitoes.
- The well-established [SkeeterBuster](https://github.com/helmingstay/SkeeterBuster) focuses on the _Aedes aegypti_ species in order to understand insecticidal control measures such as spraying [@magori2009skeeter;@Gunning2022].  A stochastic, mechanistic approach is employed.
- [OpenMalaria](https://github.com/SwissTPH/openmalaria/wiki) [@Smith2006] is an open-source C++ program enabling simulation of malaria epidemiology, typically at the village scale, in order to assess the efficacy of non-genetic malaria interventions.
- IDM's [EMOD](https://docs.idmod.org/projects/emod-malaria/en/latest/index.html) software can simulate malaria epidemiology using an agent-based approach, with spatial structure based on a network [@10.1093/femspd/fty059].  Less emphasis is spent on genetic modifications, and a single mosquito species is the focus.
- The [dynamAedes](https://cran.r-project.org/web/packages/dynamAedes/vignettes/dynamAedes_02_local.html) R package can be used to study the spatio-temporal evolution of a single mosquito species, with particular attention paid to the impact of temperature heterogeneity.
- The [exDE](https://dd-harp.github.io/exDE/) R package solves models of mosquito-borne pathogen dynamics and control [@Wu2023].  Attention is payed to sophisticated representations of mosquito lifecycles, including exogenous forcing by weather and vector control, as well as mosquito-malaria-human interatctions.  Although a single mosquito species is the focus, the framework allows for multiple species.  The code centers on traditional vector controls, rather than genetic controls.
- Berkeley's [MGDrivE](https://marshalllab.github.io/MGDrivE/) is an open-source framework to study gene-drives in mosquito populations [@mgdrive2019;@mondal2024mgdrive], which is written in R.  With regards to lifecycle dynamics, MGDrivE has similar functionality to `Mozzie`, although MGDrivE focuses on single species, in contrast to `Mozzie` where transfer of genetic constructs between (sub-)species is of interest.  MGDrivE's spatial structure is based on a network, where each node in the network could be thought of as a household, house block, or even a city.  In terms of functionality, MGDrivE is the most similar to `Mozzie`, although the numerical methods employed are quite different.

In addition to these publicly-available codes, many academic articles consider the lifecycle and spatial spread of mosquitoes, for example [@LUTAMBI2013198;@DUFOURD20131695;@north2013;@roques2016;@YAMASHITA201890;@Smith2018;@endo2018;@fc2018;@yamashita2018;@silva2020;@FANG2020149;@Bruzzone2022;@dye2024], but few have published their code.  Most appear to rely on unpublished scripts in codes such as MATLAB [@LUTAMBI2013198;@yamashita2018;@fc2018], or concentrate on specialised scenarios [@roques2016;@Bruzzone2022].

If spatially explicit, the aforementioned codes model spatial structure using a network.  In contrast, `Mozzie` uses a continuous-space (diffusion-advection equation) approach, deliberately incorporating long-range dispersal in a way that is ecologically interpretable [@Hosack2023].  In addition, `Mozzie` does not focus on single species, but concentrates on the interaction of multiple (sub-)species.  Many of the aforementioned alternatives contain human agents and the malaria parasite, which `Mozzie` does not.

# Statement of need

`Mozzie` is designed to solve problems involving:

- interacting (sub-)species of mosquitoes, with
- complicated lifecycle dynamics including transfers of genetic modifications between the (sub-)species, in
- spatially-extensive settings (such as continental scales) including the spatio-temporal dispersal of individuals (such as advection via wind).

It is anticipated that users of `Mozzie` will be researchers interested in such aspects.

Importantly, the numerical implementation of Mozzie is:

- ecologically interpretable [@Hosack2023],
- computationally and I/O efficient, and
- well tested.

This allows rapid simulation at continental scales to investigate sensitivity to input parameters, as required in risk assessments.  It is written in [Cython](https://cython.org) [@behnel2010cython] (a mixture of Python and C), and simulations are run using Python. The test coverage of the `Mozzie` codebase is over 99%, meaning it is also suitable for risk assessments that could be subject to considerable scrutiny.

# Demonstration

\autoref{demo} shows results from a `Mozzie` simulation using 2 interbreeding, hybridising and competing mosquito species [@beetonplos].  A water body runs through the domain, effectively separating it into two parts.  Mosquitoes can only traverse the water body via wind advection.  A gene-drive is introduced into one species, and the modified individuals are released at the center of the domain.  These modified individuals breed with wild mosquitoes, spreading the gene-drive. Over time, resistance forms and the genetic modification becomes ineffective in some individuals.  The wind blows in predominantly north-east or south-west directions, depending on the season, and it is assumed that a certain proportion of mosquitoes advect for 2 hours each night, which explains the clusters of rapid invasion shown in \autoref{demo}.  The code is found in the `demo` directory of the repository.

![The time taken for genetically-modified mosquito species to invade regions of the domain, after being released from the central point.  The water body is shown in cyan.\label{demo}](demo/demo.pdf)

# Acknowledgements

This work was funded by a grant to the Foundation for the National Institutes of Health from the Bill \& Melinda Gates Foundation (INV-008525).

# References
