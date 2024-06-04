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
    orcid: 0000-0000-0000-0000
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Maud El-Hachem
    orcid: 0000-0002-9312-6109
    equal-contrib: false # (This is how you can denote equal contributions between multiple authors)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Keith R Hayes
    orcid: 0000-0000-0000-0000
    equal-contrib: true # (This is how you can denote equal contributions between multiple authors)
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

`Mozzie` enables simulation of the lifecycle and spatio-temporal spread of mosquitoes.  Mosquitoes transmit malaria, and `Mozzie` is designed to support risk assessments associated with genetic alterations of mosquitoes aimed at eliminating malaria.  Treatment of the mosquito lifecycle includes competition and breeding between multiple (sub-)species, allowing researchers to assess potential transfer of the genetic alteration to other (sub-)species.  To assess the spatio-temporal spread of wild-type and genetically-altered mosquitoes, `Mozzie` can simulate the diffusion of individuals through the landscape, and long-distance dispersal via wind.  `Mozzie` can be used for risk assessments at the local, regional or continental scale.  Wind patterns, the capacity of the landscape to support mosquitoes, and so forth, can vary spatially and temporally, reflecting seasonality, local conditions and/or daily variations, etc.  The genetic alterations can include gene drives.

`Mozzie` has been used by the authors and collaborators in academic studies of a population-modifying gene-drive in malaria vectors in sub-Saharan Africa [@beetonplos] ... Geoff/Nick/Keith - please write a few words about risk-assessment work, and include refs.


# Statement of need

> JOSS says: "A Statement of need section that clearly illustrates the research purpose of the software and places it in the context of related work."

`Mozzie` is designed to be computationally and I/O efficient because risk assessments require running many simulations, varying over distributions of input parameters.  It is written in [Cython](https://cython.org) [@behnel2010cython] (a mixture of Python and C), and simulations are run using Python.  `Mozzie` is based on differential equations (or delay differential equations) that control the spatio-temporal spread and lifecycle dynamics of populations of mosquitoes.   The diffusion-advection equation, possibly including long-distance dispersal via wind, is used to model spatial spread.  Various options for lifecycle dynamics are provided, which reflect the evolving understanding of mosquito dynamics and the increasing sophistication of genetic techniques.  These include, amongst others, a logistic equation, a hybridising 2-species model [@BEETON2020110072], and a delay differential equation [@El-Hachem2024] with $N$ hybridising and competing (sub-)species, along with genetic alterations **cite**.  The test coverage of the `Mozzie` codebase is over 99%, meaning it is suitable for risk assessments of controversal measures that could be subjected to considerable scrutiny.

Alternatives to `Mozzie` include:

- MIT's [HYDREMATS](https://eltahir.mit.edu/models/hydremats/) software [@bomblies2009].  HYDREMATS is a coupled hydrology and entomology model that uses an agent-based approach for mosquito-human dynamics, and focuses on high-resolution village-scale understanding of malaria without genetically-modified mosquitoes.
- [OpenMalaria](https://github.com/SwissTPH/openmalaria/wiki) [@Smith2006] is an open-source C++ program enabling simulation of malaria epidemiology, typically at the village scale, in order to assess the effacy of non-genetic malaria interventions.
- IDM's [EMOD](https://docs.idmod.org/projects/emod-malaria/en/latest/index.html) software can simulate malaria epidemiology using an agent-based approach, with spatial structure based on a network [@10.1093/femspd/fty059].  Less emphasis is spent on genetic modifications, and a single mosquito species is the focus.
- The [dynamAedes](https://cran.r-project.org/web/packages/dynamAedes/vignettes/dynamAedes_02_local.html) R package can be used to study the spatio-temporal evolution of a single mosquito species, with particular attention paid to the impact of temperature heterogeneity.
- Berkeley's [MGDrivE](https://marshalllab.github.io/MGDrivE/) is an open-source framework to study gene-drives in mosquito populations [@mgdrive2019], which is written in R.  With regards to lifecycle dynamics, MGDrivE has similar functionality to `Mozzie`, although MGDrivE focusses on single species, in contrast to `Mozzie` where transfer of genetic constructs between (sub-)species is of prime interest.  MGDrivE's spatial structure is based on a network, where each node in the network could be thought of as a household, house block, or even a city.  In terms of functionality, MGDrivE is the most similar to `Mozzie`, although the numerical methods employed are quite different.

These model spatial structure using a network.  In contrast, `Mozzie` uses a continuous-space (diffusion-advection equation) approach, deliberately incorporating long-range dispersal.  In addition, `Mozzie` does not focus on single species, but concentrates on the interaction of multiple (sub-)species (that may or may not contain genetic modification), via competition, breeding and hybridisation.  Both long-range dispersal and gene transfer between (sub-)species are crucial in risk assessments.  Conversely, `Mozzie` does not contain human agents, nor does it consider the malaria parasite.

Many academic articles consider the lifecycle and spatial spread of mosquitoes, for example [@LUTAMBI2013198,@DUFOURD20131695,@north2013,@roques2016,@YAMASHITA201890,@Smith2018,@endo2018,@fc2018,@yamashita2018,@silva2020,@FANG2020149,@Bruzzone2022,@dye2024], but few have published their code.  Most appear to rely on unpublished scripts in codes such as MATLAB [@LUTAMBI2013198,@yamashita2018,@fc2018], or concentrate on specialised scenarios [@roques2016,@Bruzzone2022].

** NICE FIGURE **

# Other sections and info

- A list of key references, including to other software addressing related needs. Note that the references should include full names of venues, e.g., journals and conferences, not abbreviations only understood in the context of a specific discipline.
- Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.
- Acknowledgement of any financial support.

# Acknowledgements

This work was supported, in whole or in part, by the Bill & Melinda Gates Foundation (**grant number**).  Under the grant conditions of the Foundation, a Creative Commons Attribution 4.0 Generic License has already been assigned to the Author Accepted Manuscript version that might arise from this submission.

See [instructions](https://joss.readthedocs.io/en/latest/paper.html) for instructions about formulae, references, tables, figures, markdown, footnotes, headings, etc.  See [example](https://joss.readthedocs.io/en/latest/example_paper.html) for an example.

# References
  