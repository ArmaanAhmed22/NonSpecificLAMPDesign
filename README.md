# Sensitive LAMP Primer Design
Loop-mediated isothermal AMPlification has been used in conjunction with Cas12b for point-of-care detection (called POC SHERLOCK). Here, I design a pipeline to create sensitive LAMP primers for POC SHERLOCK, specifically against HIV-1 protease.

### Algorithm
1. Get all possible primer sets
2. Score and filter out primer sets based on primer length, Tm, dG, and spacing
3. With the remaining, efficacious primer sets, keep the ones with the highest conservativety
4. Repeat algorithm to create primer sets that targets what is not currently targetable (repeat until overall sensitivity above a certain threshold)