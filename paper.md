---
title: 'Mozzie: a computationally efficient simulator for the spatio-temporal modelling of mosquitoes'
tags:
  - Python
  - Cython
  - malaria
  - mosquito
  - risk analysis
authors:
  - name: Andy Wilkins
    orcid: 0000-0001-6472-9560
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Nicholas J. Beeton
    orcid: 0000-0000-0000-0000 # TODO
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Maud El-Hachem
    orcid: 0000-0002-9312-6109
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Keith R Hayes
    orcid: 0000-0003-4094-3575
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Geoffrey R. Hosack
    orcid: 0000-0002-6462-6817
    equal-contrib: false
    affiliation: "1"
affiliations:
 - name: Commonwealth Scientific and Industrial Organisation (CSIRO), Australia
   index: 1
date: 1 June 2024
bibliography: paper.bib
---

The paper should be between 250-1000 words. Authors submitting papers significantly longer than 1000 words may be asked to reduce the length of their paper.

# Summary

> JOSS says: "A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.  We also require that authors explain the research applications of the software."

`Mozzie` enables simulation of the lifecycle and spatial spread of mosquitoes.  `Mozzie` can be used to assess risks associated with disease-control strategies at local, regional or continental scales, most particularly, strategies involving genetic alterations of mosquitoes to eliminate malaria.

`Mozzie` simulates a population-dynamics model that uses differential equations or delay differential equations [@BOHNER2018114;@El-Hachem2024] to describe the spread and persistence of genetic alterations that modify or suppress populations of disease-carrying mosquitoes.  Genetic alterations are flexibly modelled: these can involve any number of alleles; Mendelian or non-Mendelian inheritance, including gene drives; they can be self-limiting or self-sustaining, and can include the emergence of resistant allelles.  The model allows simulation of $N$ mosquito species.  It incorporates mate-choice, hybridisation and intra-specific competition that occur within complexes of sibling species [@BEETON2020110072], such as the three main malaria vectors in the _Anopheles gambiae_ complex.  This fills a gap that currently exists among similar models, allowing researchers to assess potential transfer of the genetic alterations between (sub-)species.  `Mozzie` supports age- or stage-specific population structure, spatial and temporal variations in birth rate, local diffusion and wind-assisted, long range, advection.  Spatio-temporal heterogeneity is supported.  For example, wind patterns and the capacity of the landscape to support mosquitoes can vary spatially and temporally, reflecting seasonality, local conditions and/or daily variations, etc.  

Conversely, `Mozzie` does not contain human agents, nor does it consider the effect of genetic control strategies on the prevalence of pathogens such as the malaria parasite, among human or animal populations.

`Mozzie` has been used by the authors to simulate the spread across sub-Saharan Africa of a theoretical, population-modifying, gene drive in _Anopheles gambiae s.s._ and _Anopheles coluzzii_, and to predict the spread of Target Malaria's Paternal Male Bias construct [@Galizi2014] following a proposed field-release of genetically modified _Anopheles coluzzi_ male mosquitoes in Burkina Faso [@Hosack2023].


AN ALTERATE PARAGRAPH:

`Mozzie` enables simulation of the lifecycle and spatial spread of mosquitoes.  Mosquitoes transmit malaria, and `Mozzie` is designed to support risk assessments associated with genetic alterations of mosquitoes aimed at eliminating malaria.  Treatment of the mosquito lifecycle includes competition and breeding between multiple (sub-)species, allowing researchers to assess potential transfer of the genetic alteration to other (sub-)species.  To assess the spatial spread of wild-type and genetically-altered mosquitoes, `Mozzie` can simulate the local diffusion of individuals through the landscape, and long-distance dispersal via wind.  `Mozzie` can be used for risk assessments at the local, regional or continental scale.  For example, wind patterns and the capacity of the landscape to support mosquitoes can vary spatially and temporally, reflecting seasonality, local conditions and/or daily variations, etc.  Genetic alterations can be flexibly modelled with any number of alleles and the capacity for Mendelian or non-Mendelian inheritance, including gene drives.

# Statement of need

> JOSS says: "A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work."

The numerical implementation of Mozzie is ecologically interpretable [@Hosack2023] and computationally and I/O efficient.  This allows simulations at continental scales and the investigation of varying input parameters, as required in risk assessments.  It is written in [Cython](https://cython.org) [@behnel2010cython] (a mixture of Python and C), and simulations are run using Python. The test coverage of the `Mozzie` codebase is over 99%, meaning it is also suitable for risk assessments of vector control strategies that could be subject to considerable scrutiny.

Alternatives to `Mozzie` include:

- MIT's [HYDREMATS](https://eltahir.mit.edu/models/hydremats/) software [@bomblies2009].  HYDREMATS is a coupled hydrology and entomology model that uses an agent-based approach for mosquito-human dynamics, and focuses on high-resolution village-scale understanding of malaria without genetically-modified mosquitoes.
- The well-established [SkeeterBuster](https://github.com/helmingstay/SkeeterBuster) focuses on the _Aedes Aegypti_ species in order to understand insecticidal control measures such as spraying [@magori2009skeeter;@Gunning2022].  A stochastic, mechanistic approach is employed.
- [OpenMalaria](https://github.com/SwissTPH/openmalaria/wiki) [@Smith2006] is an open-source C++ program enabling simulation of malaria epidemiology, typically at the village scale, in order to assess the efficacy of non-genetic malaria interventions.
- IDM's [EMOD](https://docs.idmod.org/projects/emod-malaria/en/latest/index.html) software can simulate malaria epidemiology using an agent-based approach, with spatial structure based on a network [@10.1093/femspd/fty059].  Less emphasis is spent on genetic modifications, and a single mosquito species is the focus.
- The [dynamAedes](https://cran.r-project.org/web/packages/dynamAedes/vignettes/dynamAedes_02_local.html) R package can be used to study the spatio-temporal evolution of a single mosquito species, with particular attention paid to the impact of temperature heterogeneity.
- The [exDE](https://dd-harp.github.io/exDE/) R package solves models of mosquito-borne pathogen dynamics and control [@Wu2023].  Attention is payed to sophisticated representations of mosquito lifecycles, including exogenous forcing by weather and vector control, as well as mosquito-malaria-human interatctions.  Although a single mosquito species is the focus, the framework allows for multiple species.  The code centers on traditional vector controls, rather than genetic controls.
- Berkeley's [MGDrivE](https://marshalllab.github.io/MGDrivE/) is an open-source framework to study gene-drives in mosquito populations [@mgdrive2019;@mondal2024mgdrive], which is written in R.  With regards to lifecycle dynamics, MGDrivE has similar functionality to `Mozzie`, although MGDrivE focuses on single species, in contrast to `Mozzie` where transfer of genetic constructs between (sub-)species is of prime interest.  MGDrivE's spatial structure is based on a network, where each node in the network could be thought of as a household, house block, or even a city.  In terms of functionality, MGDrivE is the most similar to `Mozzie`, although the numerical methods employed are quite different.

If spatially explicit, these model spatial structure using a network.  In contrast, `Mozzie` uses a continuous-space (diffusion-advection equation) approach, deliberately incorporating long-range dispersal in a way that is ecologically interpretable [@Hosack2023].  In addition, `Mozzie` does not focus on single species, but concentrates on the interaction of multiple (sub-)species (that may or may not contain genetic modification).  Many of the aforementioned alternatives contain human agents and the malaria parasite, which `Mozzie` does not.

In addition to the publicly available codes mentioned above, many academic articles consider the lifecycle and spatial spread of mosquitoes, for example [@LUTAMBI2013198;@DUFOURD20131695;@north2013;@roques2016;@YAMASHITA201890;@Smith2018;@endo2018;@fc2018;@yamashita2018;@silva2020;@FANG2020149;@Bruzzone2022;@dye2024], but few have published their code.  Most appear to rely on unpublished scripts in codes such as MATLAB [@LUTAMBI2013198,@yamashita2018,@fc2018], or concentrate on specialised scenarios [@roques2016;@Bruzzone2022].

** NICE FIGURE **
# Maud's suggestion: maybe a UML Component Diagram, there is an nice example in this paper https://joss.theoj.org/papers/10.21105/joss.06384
# Andy was thinking of a figure from our risk assessment work
TODO

# Other sections and info

- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.
- Acknowledgement of any financial support.

# Acknowledgements

This work was supported, in whole or in part, by the Bill & Melinda Gates Foundation (**grant number**) TODO.  Under the grant conditions of the Foundation, a Creative Commons Attribution 4.0 Generic License has already been assigned to the Author Accepted Manuscript version that might arise from this submission.

See [instructions](https://joss.readthedocs.io/en/latest/paper.html) for instructions about formulae, references, tables, figures, markdown, footnotes, headings, etc.  See [example](https://joss.readthedocs.io/en/latest/example_paper.html) for an example.

# References
