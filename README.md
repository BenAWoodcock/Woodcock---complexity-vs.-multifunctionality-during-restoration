# Woodcock---complexity-vs.-multifunctionality-during-restoration
This analysis assesses the relationship between ecological complexity and system multifunctionality during habitat restoration of grasslands and woodlands.

**Code and data for analyses presented in ‘Restoration of ecological complexity enhances system multifunctionality’ Woodcock et al**

**Abstract**
Habitat loss is a major driver of biodiversity decline that has contributed to the degradation of key ecosystem functions.  Restoration of lost and degraded habitats has been used to reverse these trends.  While a focus on replicating reference communities has been widely used in the past, this emphasis on the re-establishment of expected species may limit our ability to effectively restore critical functions under adverse environmental drivers like climate change.  Restoration of ecological complexity, which considers the number of system components, connections and feedback amongst them, provides an alternative endpoint that has a stronger theoretical expectation to support higher functioning ecosystems. Here we assess this hypothesis by deriving aggregate measures of multifunctionality across 18 grassland and 18 woodland restoration sites in the UK.  Multifunctionality is based on measures of soil functions, carbon stocks, net primary productivity, decomposition, herbivory, pollination, and predation.   We also derive measures of ecological complexity for each site based on bacterial, fungal, plant and insect communities.  We show consistent and clear relationships between ecological complexity and multifunctionality across grasslands and woodlands, both individually and when combined.  In contrast, plant species richness represents a poor predictor of ecological multifunctionality.  This emphasises the need for new metrics to support management decisions during restoration if multifunctionality is to be preserved under changing environmental pressures.  We propose that the remit of restoration, and with it the basis for quantifying success, needs to be expanded to consider wider complexity endpoints where high functioning, resilient systems are desired.  

Below I list files contained in this repository combined with meta-data on individual variables.

**RestREco_KP2_Complex_vs_Function.R**
R statistical package data file containing code for analyses presented in paper.

**Site_info.csv**
Habitat =  grassland or woodland
Region = England (grassland) or Scotland (woodland)
FLU =  former land use
Age =  age of restoration in years
Latitude =  site latitude
Longitude =  site longitude

**Compelxity_metrics.csv**
Site.code = unique code for each of the 18 woodland and 18 grassland sites
DD.early.mean =  Grassland mean sward height taken early in the spring	
DD.late.mean	=  Grassland mean sward height taken early in the summer	
trees.DBH.mean = woodland tree mean diameter at breast height
trees.DBH.sd = woodland tree standard deviation in diameter at breast height
trees.stem.dens.invert = density of trees	
trees.BA.ha	= Basal area of trees per hectare
Plant.sp.incl.sp =  plant species richness resolved to the species level
Plant.count.incl.sp	=  Site-level relative abundance of plants calculated as the sum species occurrences within the ten quadrats
Plant.FD.rich	= plant community functional richness
Plant.FD.div	= plant community functional diversity
Invert.Order.N = total invertebrate abundance	
Invert.sp.rich	= total invertebrate species richness
Invert.Biomass =  invertebrate biomass
Invert.Herb.N	= abundance of invertebrate herbivores
Invert.Omni.N= abundance of invertebrate omnivores	
Invert.Det.N= abundance of invertebrate detritivores	
Invert.Pred.N= abundance of invertebrate predators	
Invert.FD.rich.sp =  invertebrate functional richness	
Invert.FD.div.sp =  invertebrate functional diversity	
connectance	= bipartite web connectance for plant-pollinator (grassland) and plant-herbivore (woodland) invertebrate communities
Generality= bipartite web generality for plant-pollinator (grassland) and plant-herbivore (woodland) invertebrate communities	
Vulnerability= bipartite web vulnerability for plant-pollinator (grassland) and plant-herbivore (woodland) invertebrate communities	
Bacteria.H	 = soil  bacterial diversity
Fungi.H =  soil fungal diversity

**Process_metrics.csv**
Site.code = unique code for each of the 18 woodland and 18 grassland sites
Soil_resp_ave	 =  Summed carbon dioxide emissions from bacterial processing of an organic substrate.
Soil_resp_SD	= Standard deviation in carbon dioxide emissions from bacterial processing of an organic substrate.
Fungi.FD.Chao1 = Soil fungal Chao1 diversity of functional gene sequences involved in processing organic and inorganic substrates in soils.  	
Bacteria.FD.Chao1	= Soil bacterial Chao1 diversity of functional gene sequences involved in processing organic and inorganic substrates in soils.  	
Soil.total.organic.C =  Total soil organic carbon (TOC) determined by combustion following removal of inorganic soil carbonates.
BAIT.LAMINA.high.perc.removal.Horisonal =  percentage removal of carbon substrate from bait-lamina strips used to assess the loss of an coarse organic substrate.  Soil surface assessment.	
BAIT.LAMINA.high.perc.removal.Vertical =  =  percentage removal of carbon substrate from bait-lamina strips used to assess the loss of an coarse organic substrate.  Soil sub-surface assessment.	
Biomass.kg.ha =  Above ground biomass accumulation in grassland swards (within year) and tree biomass (accumulated across multiple years) for woodlands.	
eDNA_diet_Fam_PtMad =  Family richness of invertebrate prey fed upon by the ground beetle Pterositcus madidus as recorded form all sites.	
Ave.Nbites.small.beetles	= Attack rates in terms of average number of bites from small beetles on model plasticine caterpillars.
Ave.Nbites.large.beetles	= Attack rates in terms of average number of bites from large beetles on model plasticine caterpillars.
Ave.Nbites.rodent = Attack rates in terms of average number of bites from rodents on model plasticine caterpillars.
Ave.Nbites.bird = Attack rates in terms of average number of bites from small beetles on model plasticine caterpillars.
Ave.Nbites.total.bites	= Attack rates in terms of average number of bites from all taxa on model plasticine caterpillars.
Prop.bitten.small.beetles	= Attack rates in terms of the proportion of plasticine caterpillars that were attacked by small beetles.
Prop.bitten.large.beetles	= Attack rates in terms of the proportion of plasticine caterpillars that were attacked by large beetles.
Prop.bitten.rodent	= Attack rates in terms of the proportion of plasticine caterpillars that were attacked by rodents.
Prop.bitten.bird	= Attack rates in terms of the proportion of plasticine caterpillars that were attacked by birds.
Poll.Leontodon.seed.set	 =  Seed set of Leontodon flowers in grasslands attributable to insect pollinators.
Poll.wood.seed.set		 =  Seed set of Rannunculus flowers in woodlands attributable to insect pollinators.
Perc.herb.Oak =  Rates of herbivory on leaves of oak in woodlands.
Perc.herb.Birch =  Rates of herbivory on leaves of birch in woodlands.	
Perc.herb.Plantago=  Rates of herbivory on leaves of Plantago in grasslands.
Perc.herb.Ran =  Rates of herbivory on leaves of Rannunculus in grasslands.
