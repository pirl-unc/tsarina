# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""
Geographic HLA allele frequency data for region-weighted vaccine design.

Provides curated per-region HLA class I allele frequencies drawn from
population-level studies and bone marrow registries, along with 2024
estimated regional populations.  These data support region-weighted
coverage calculations when selecting peptide-vaccine allele panels.
"""

from __future__ import annotations

from collections.abc import Iterable

import pandas as pd

REGION_PRIORITY_ROWS: list[dict] = [
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "A",
        "allele": "HLA-A*02:01",
        "frequency": 0.2280,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "A",
        "allele": "HLA-A*24:02",
        "frequency": 0.1230,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "A",
        "allele": "HLA-A*01:01",
        "frequency": 0.1150,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "B",
        "allele": "HLA-B*51:01",
        "frequency": 0.0980,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "B",
        "allele": "HLA-B*18:01",
        "frequency": 0.0950,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "B",
        "allele": "HLA-B*35:01",
        "frequency": 0.0800,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Representative B*35 family allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "C",
        "allele": "HLA-C*04:01",
        "frequency": 0.1750,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "C",
        "allele": "HLA-C*07:01",
        "frequency": 0.1710,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "Europe",
        "proxy": "Italy donor registry",
        "locus": "C",
        "allele": "HLA-C*06:02",
        "frequency": 0.0990,
        "resolution": "exact",
        "source_label": "Italian Bone Marrow Donor Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC6771744/",
        "note": "Top European proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "A",
        "allele": "HLA-A*02:01",
        "frequency": 0.1589,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "A",
        "allele": "HLA-A*24:02",
        "frequency": 0.0937,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "A",
        "allele": "HLA-A*01:01",
        "frequency": 0.0773,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "A",
        "allele": "HLA-A*11:01",
        "frequency": 0.0648,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "A",
        "allele": "HLA-A*68:01",
        "frequency": 0.0579,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "B",
        "allele": "HLA-B*50:01",
        "frequency": 0.1202,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab-specific B allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "B",
        "allele": "HLA-B*51:01",
        "frequency": 0.1049,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "B",
        "allele": "HLA-B*08:01",
        "frequency": 0.0723,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "B",
        "allele": "HLA-B*40:06",
        "frequency": 0.0361,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Important South Asia / Arab overlap allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "B",
        "allele": "HLA-B*58:01",
        "frequency": 0.0319,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*06:02",
        "frequency": 0.1396,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*04:01",
        "frequency": 0.1320,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1289,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Top Arab proxy allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*15:02",
        "frequency": 0.0920,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*07:01",
        "frequency": 0.0729,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*17:01",
        "frequency": 0.0538,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "MENA allele worth considering in an extended panel.",
    },
    {
        "region": "MENA / Arab",
        "proxy": "Kuwait",
        "locus": "C",
        "allele": "HLA-C*12:03",
        "frequency": 0.0513,
        "resolution": "exact",
        "source_label": "Kuwaiti population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7024892/",
        "note": "Common Arab allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "A",
        "allele": "HLA-A*11:01",
        "frequency": 0.2680,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Top East Asian A allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "A",
        "allele": "HLA-A*24:02",
        "frequency": 0.1900,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Top East Asian A allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "A",
        "allele": "HLA-A*02:07",
        "frequency": 0.0940,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "East Asia-enriched A allele commonly missed in European-centric panels.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "A",
        "allele": "HLA-A*02:01",
        "frequency": 0.0860,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Common East Asian A allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "A",
        "allele": "HLA-A*02:06",
        "frequency": 0.0520,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Common East Asian A allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "A",
        "allele": "HLA-A*02:03",
        "frequency": 0.0400,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Common East Asian A allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "B",
        "allele": "HLA-B*46:01",
        "frequency": 0.1440,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Major East Asian B allele excluded from most European-focused work.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "B",
        "allele": "HLA-B*40:01",
        "frequency": 0.1170,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Top East Asian B allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "B",
        "allele": "HLA-B*13:01",
        "frequency": 0.0840,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Common East Asian B allele not in most compact panels.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "C",
        "allele": "HLA-C*01:02",
        "frequency": 0.2030,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Major East Asian C allele excluded from most European-focused panels.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1530,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Common East Asian C allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Yunnan Han",
        "locus": "C",
        "allele": "HLA-C*03:04",
        "frequency": 0.1170,
        "resolution": "exact",
        "source_label": "Yunnan Han population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4068047/",
        "note": "Common East Asian C allele.",
    },
    {
        "region": "East Asia",
        "proxy": "Singapore Chinese",
        "locus": "A",
        "allele": "HLA-A*11:01",
        "frequency": 0.2827,
        "resolution": "exact",
        "source_label": "Singapore Chinese donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent East Asian validation of A*11:01 dominance.",
    },
    {
        "region": "East Asia",
        "proxy": "Singapore Chinese",
        "locus": "B",
        "allele": "HLA-B*40:01",
        "frequency": 0.1916,
        "resolution": "exact",
        "source_label": "Singapore Chinese donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent East Asian validation of B*40:01.",
    },
    {
        "region": "East Asia",
        "proxy": "Singapore Chinese",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1996,
        "resolution": "exact",
        "source_label": "Singapore Chinese donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent East Asian validation of C*07:02.",
    },
    {
        "region": "East Asia",
        "proxy": "China South Han pop 2",
        "locus": "C",
        "allele": "HLA-C*04:03",
        "frequency": 0.0132,
        "resolution": "exact",
        "source_label": "Allele Frequency Net Database: China South Han pop 2",
        "source_url": "https://www.allelefrequencies.net/hla6006b.asp?hla_population=2768",
        "note": "Specific-population frequency support for a Sarkizova frequent HLA-C allotype.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "A",
        "allele": "HLA-A*11:01",
        "frequency": 0.2606,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Top Southeast Asian A allele.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "B",
        "allele": "HLA-B*46:01",
        "frequency": 0.1404,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Top Southeast Asian B allele.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "B",
        "allele": "HLA-B*15:02",
        "frequency": 0.0766,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Major Southeast Asian B allele often missed in compact global panels.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "B",
        "allele": "HLA-B*40:01",
        "frequency": 0.0660,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Common Southeast Asian B allele.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "B",
        "allele": "HLA-B*58:01",
        "frequency": 0.0638,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Common Southeast Asian B allele.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "C",
        "allele": "HLA-C*01:02",
        "frequency": 0.1713,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Top Southeast Asian C allele.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1191,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Common Southeast Asian C allele.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Thailand",
        "locus": "C",
        "allele": "HLA-C*08:01",
        "frequency": 0.1032,
        "resolution": "exact",
        "source_label": "Thai population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC7057685/",
        "note": "Important Southeast Asian C allele not in compact baseline panels.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Singapore Malay",
        "locus": "A",
        "allele": "HLA-A*11:01",
        "frequency": 0.1804,
        "resolution": "exact",
        "source_label": "Singapore Malay donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent Southeast Asian validation of A*11:01.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Singapore Malay",
        "locus": "B",
        "allele": "HLA-B*15:02",
        "frequency": 0.1148,
        "resolution": "exact",
        "source_label": "Singapore Malay donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent Southeast Asian validation of B*15:02.",
    },
    {
        "region": "Southeast Asia",
        "proxy": "Singapore Malay",
        "locus": "C",
        "allele": "HLA-C*08:01",
        "frequency": 0.2008,
        "resolution": "exact",
        "source_label": "Singapore Malay donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent Southeast Asian validation of C*08:01.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "A",
        "allele": "HLA-A*01:01",
        "frequency": 0.2541,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for dominant HLA-A*01 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "A",
        "allele": "HLA-A*02:01",
        "frequency": 0.2483,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for dominant HLA-A*02 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "A",
        "allele": "HLA-A*11:01",
        "frequency": 0.1753,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for dominant HLA-A*11 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "A",
        "allele": "HLA-A*24:02",
        "frequency": 0.1027,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-A*24 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "A",
        "allele": "HLA-A*03:01",
        "frequency": 0.0907,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-A*03 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "B",
        "allele": "HLA-B*35:01",
        "frequency": 0.2054,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-B*35 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "B",
        "allele": "HLA-B*15:01",
        "frequency": 0.1536,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-B*15 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "B",
        "allele": "HLA-B*40:01",
        "frequency": 0.1359,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-B*40 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "B",
        "allele": "HLA-B*07:02",
        "frequency": 0.1014,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-B*07 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "B",
        "allele": "HLA-B*44:02",
        "frequency": 0.0779,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-B*44 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.2806,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-C*07 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "C",
        "allele": "HLA-C*04:01",
        "frequency": 0.2042,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-C*04 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "C",
        "allele": "HLA-C*03:04",
        "frequency": 0.1555,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-C*03 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "C",
        "allele": "HLA-C*06:02",
        "frequency": 0.1304,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-C*06 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "North India",
        "locus": "C",
        "allele": "HLA-C*12:03",
        "frequency": 0.0527,
        "resolution": "proxy_group",
        "source_label": "North Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC12392222/",
        "note": "Representative high-resolution allele for HLA-C*12 antigen group.",
    },
    {
        "region": "South Asia",
        "proxy": "Singapore Indians",
        "locus": "A",
        "allele": "HLA-A*24:02",
        "frequency": 0.1659,
        "resolution": "exact",
        "source_label": "Singapore Indian donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent South Asian proxy allele.",
    },
    {
        "region": "South Asia",
        "proxy": "Singapore Indians",
        "locus": "B",
        "allele": "HLA-B*40:06",
        "frequency": 0.1008,
        "resolution": "exact",
        "source_label": "Singapore Indian donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "South Asia-enriched B allele usually absent from European-centric panels.",
    },
    {
        "region": "South Asia",
        "proxy": "Singapore Indians",
        "locus": "C",
        "allele": "HLA-C*06:02",
        "frequency": 0.1396,
        "resolution": "exact",
        "source_label": "Singapore Indian donors",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9873421/",
        "note": "Independent South Asian proxy allele.",
    },
    {
        "region": "South Asia",
        "proxy": "17th IHIW Indian cohort",
        "locus": "A",
        "allele": "HLA-A*01:01",
        "frequency": 0.1650,
        "resolution": "exact",
        "source_label": "17th IHIW Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8315142/",
        "note": "India-specific A allele frequency reported in global workshop summary.",
    },
    {
        "region": "South Asia",
        "proxy": "17th IHIW Indian cohort",
        "locus": "A",
        "allele": "HLA-A*02:11",
        "frequency": 0.0710,
        "resolution": "exact",
        "source_label": "17th IHIW Indian population",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8315142/",
        "note": "South Asia-specific A allele elevated in Asian Indians.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "A",
        "allele": "HLA-A*30:01",
        "frequency": 0.1310,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-A*30 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "A",
        "allele": "HLA-A*68:02",
        "frequency": 0.1280,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-A*68 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "A",
        "allele": "HLA-A*02:01",
        "frequency": 0.1230,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-A*02 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "B",
        "allele": "HLA-B*15:03",
        "frequency": 0.1260,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative SSA-focused allele for dominant HLA-B*15 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "B",
        "allele": "HLA-B*58:02",
        "frequency": 0.1170,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative SSA-focused allele for dominant HLA-B*58 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "B",
        "allele": "HLA-B*44:03",
        "frequency": 0.0990,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-B*44 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1680,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-C*07 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "C",
        "allele": "HLA-C*06:02",
        "frequency": 0.1640,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-C*06 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Free State African descent",
        "locus": "C",
        "allele": "HLA-C*04:01",
        "frequency": 0.1280,
        "resolution": "proxy_group",
        "source_label": "Free State African descent group",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC8085688/",
        "note": "Representative high-resolution allele for HLA-C*04 group.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "South African Bone Marrow Registry",
        "locus": "B",
        "allele": "HLA-B*07:02",
        "frequency": 0.0820,
        "resolution": "exact",
        "source_label": "South African Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5937380/",
        "note": "Exact SSA registry allele.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "South African Bone Marrow Registry",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1800,
        "resolution": "exact",
        "source_label": "South African Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC5937380/",
        "note": "Exact SSA registry allele.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Southern African review",
        "locus": "A",
        "allele": "HLA-A*74:01",
        "frequency": None,
        "resolution": "qualitative",
        "source_label": "Southern African perspective review",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4549606/",
        "note": "African-enriched allele highlighted as underrepresented outside SSA studies.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Southern African review",
        "locus": "B",
        "allele": "HLA-B*57:03",
        "frequency": None,
        "resolution": "qualitative",
        "source_label": "Southern African perspective review",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4549606/",
        "note": "African-enriched allele commonly absent from compact global panels.",
    },
    {
        "region": "Sub-Saharan Africa",
        "proxy": "Southern African review",
        "locus": "B",
        "allele": "HLA-B*81:01",
        "frequency": None,
        "resolution": "qualitative",
        "source_label": "Southern African perspective review",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC4549606/",
        "note": "African-enriched allele commonly absent from compact global panels.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "A",
        "allele": "HLA-A*24:02",
        "frequency": 0.2087,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Top Latin American proxy allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "A",
        "allele": "HLA-A*02:01",
        "frequency": 0.1611,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American proxy allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "A",
        "allele": "HLA-A*01:01",
        "frequency": 0.0706,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American proxy allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "A",
        "allele": "HLA-A*03:01",
        "frequency": 0.0683,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American proxy allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "A",
        "allele": "HLA-A*29:02",
        "frequency": 0.0522,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Latin America-relevant allele worth adding in extended panels.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "B",
        "allele": "HLA-B*35:01",
        "frequency": 0.0769,
        "resolution": "proxy_group",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Representative B*35 family allele; top exact Colombian allele is B*35:43.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "B",
        "allele": "HLA-B*40:02",
        "frequency": 0.0718,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Latin America-relevant B allele worth adding in extended panels.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "B",
        "allele": "HLA-B*44:03",
        "frequency": 0.0607,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American B allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "B",
        "allele": "HLA-B*07:02",
        "frequency": 0.0536,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American B allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "C",
        "allele": "HLA-C*04:01",
        "frequency": 0.1540,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Top Latin American C allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "C",
        "allele": "HLA-C*01:02",
        "frequency": 0.1049,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American C allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "C",
        "allele": "HLA-C*07:02",
        "frequency": 0.1044,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American C allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "C",
        "allele": "HLA-C*07:01",
        "frequency": 0.0893,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American C allele.",
    },
    {
        "region": "Latin America",
        "proxy": "Colombia 2022 registry",
        "locus": "C",
        "allele": "HLA-C*03:04",
        "frequency": 0.0805,
        "resolution": "exact",
        "source_label": "Colombian Bone Marrow Registry",
        "source_url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC9869256/",
        "note": "Common Latin American C allele.",
    },
]

GLOBAL_ALLELE_FREQUENCY_ROWS: list[dict] = [
    {
        "allele": "HLA-A*23:01",
        "frequency": 0.02096,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-A*26:01",
        "frequency": 0.03352,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-A*30:02",
        "frequency": 0.00838,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-A*31:01",
        "frequency": 0.02462,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-A*32:01",
        "frequency": 0.03200,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-A*33:01",
        "frequency": 0.00749,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-B*27:05",
        "frequency": 0.02952,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-B*53:01",
        "frequency": 0.00696,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-B*57:01",
        "frequency": 0.03134,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*02:02",
        "frequency": 0.04281,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*03:02",
        "frequency": 0.00779,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*03:03",
        "frequency": 0.04237,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*05:01",
        "frequency": 0.05597,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*07:04",
        "frequency": 0.01709,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*08:02",
        "frequency": 0.02483,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*12:02",
        "frequency": 0.01526,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*14:02",
        "frequency": 0.01383,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
    {
        "allele": "HLA-C*16:01",
        "frequency": 0.02241,
        "source_label": "CIWD global allele frequency",
        "source_url": "https://www.mdpi.com/2073-4409/10/11/3048",
        "note": "Global fallback for panel coverage when no region proxy row is available.",
    },
]


_GLOBAL_ALLELE_FREQUENCY_BASE_ROWS = GLOBAL_ALLELE_FREQUENCY_ROWS


def _allele_locus(allele: str) -> str:
    return allele.removeprefix("HLA-").split("*", 1)[0]


def _global_frequency_row(
    allele: str,
    frequency: float,
    *,
    source_label: str = "CIWD global allele frequency",
    source_url: str = "https://www.mdpi.com/2073-4409/10/11/3048",
    proxy: str = "CIWD global model population",
    note: str = "Published global average for panel-frequency audit and fallback coverage.",
) -> dict:
    return {
        "region": "Global",
        "proxy": proxy,
        "locus": _allele_locus(allele),
        "allele": allele,
        "frequency": frequency,
        "resolution": "exact",
        "source_label": source_label,
        "source_url": source_url,
        "note": note,
    }


def _compatible_global_row(row: dict) -> dict:
    allele = str(row["allele"])
    return {
        "region": row.get("region", "Global"),
        "proxy": row.get("proxy", "CIWD global model population"),
        "locus": row.get("locus", _allele_locus(allele)),
        "allele": allele,
        "frequency": row["frequency"],
        "resolution": row.get("resolution", "exact"),
        "source_label": row["source_label"],
        "source_url": row["source_url"],
        "note": row["note"],
    }


_GLOBAL_CIWD_ADDITIONAL_ALLELE_FREQUENCIES: dict[str, float] = {
    "HLA-A*01:01": 0.05217,
    "HLA-A*02:01": 0.13231,
    "HLA-A*02:03": 0.00695,
    "HLA-A*02:06": 0.00901,
    "HLA-A*03:01": 0.06017,
    "HLA-A*11:01": 0.04706,
    "HLA-A*24:02": 0.08653,
    "HLA-A*29:02": 0.01041,
    "HLA-A*30:01": 0.00924,
    "HLA-A*68:01": 0.02606,
    "HLA-A*68:02": 0.01048,
    "HLA-B*07:02": 0.05529,
    "HLA-B*08:01": 0.03277,
    "HLA-B*15:01": 0.02271,
    "HLA-B*15:02": 0.00261,
    "HLA-B*18:01": 0.02336,
    "HLA-B*35:01": 0.03821,
    "HLA-B*40:01": 0.01916,
    "HLA-B*40:02": 0.01648,
    "HLA-B*44:02": 0.03374,
    "HLA-B*44:03": 0.03113,
    "HLA-B*46:01": 0.01125,
    "HLA-B*51:01": 0.05382,
    "HLA-B*58:01": 0.01254,
    "HLA-C*01:02": 0.03530,
    "HLA-C*03:04": 0.01919,
    "HLA-C*04:01": 0.06437,
    "HLA-C*06:02": 0.05748,
    "HLA-C*07:01": 0.12033,
    "HLA-C*07:02": 0.07478,
    "HLA-C*08:01": 0.02085,
    "HLA-C*12:03": 0.03936,
    "HLA-C*15:02": 0.02678,
    "HLA-C*17:01": 0.00752,
}

GLOBAL_ALLELE_FREQUENCY_ROWS = [
    _compatible_global_row(row)
    for row in [
        *_GLOBAL_ALLELE_FREQUENCY_BASE_ROWS,
        *[
            _global_frequency_row(allele, frequency)
            for allele, frequency in _GLOBAL_CIWD_ADDITIONAL_ALLELE_FREQUENCIES.items()
        ],
        _global_frequency_row(
            "HLA-C*04:03",
            0.019,
            source_label="Sarkizova HLA-C global allele frequency",
            source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12738900/",
            proxy="Sarkizova global HLA-C allotype set",
            note=(
                "Published HLA-C allotype frequency; used because CIWD Table A2 "
                "does not include HLA-C*04:03."
            ),
        ),
        _global_frequency_row(
            "HLA-C*14:03",
            0.015,
            source_label="Sarkizova HLA-C global allele frequency",
            source_url="https://pmc.ncbi.nlm.nih.gov/articles/PMC12738900/",
            proxy="Sarkizova global HLA-C allotype set",
            note="Published HLA-C allotype frequency for the legacy global51_abc panel.",
        ),
    ]
]


REGION_POPULATIONS: dict[str, float] = {
    "Europe": 743.9,
    "MENA / Arab": 599.9,
    "East Asia": 1_603.2,
    "Southeast Asia": 704.8,
    "South Asia": 2_076.8,
    "Latin America": 672.1,
    "Sub-Saharan Africa": 1_311.5,
}


def region_allele_frequencies() -> pd.DataFrame:
    """Return a DataFrame of per-region HLA allele frequencies.

    Each row corresponds to one entry in :data:`REGION_PRIORITY_ROWS`.
    """
    return pd.DataFrame(REGION_PRIORITY_ROWS)


def global_allele_frequencies() -> pd.DataFrame:
    """Return published global HLA allele-frequency averages.

    These rows use the same schema as :func:`region_allele_frequencies` with
    ``region == "Global"``. Coverage calculations prefer numeric regional
    proxies and fall back to these rows only when no regional proxy exists.
    """
    return pd.DataFrame(GLOBAL_ALLELE_FREQUENCY_ROWS)


def allele_frequency_rows(alleles: Iterable[str] | None = None) -> pd.DataFrame:
    """Return regional proxy and global-average allele frequencies together.

    Rows are on the same 0-1 allele-frequency scale. ``frequency_scope``
    distinguishes sub-population proxy rows from published global averages.
    """
    regional = region_allele_frequencies().copy()
    regional["frequency_scope"] = "regional_proxy"
    global_df = global_allele_frequencies().copy()
    global_df["frequency_scope"] = "published_global"
    out = pd.concat([regional, global_df], ignore_index=True, sort=False)
    columns = [
        "frequency_scope",
        "region",
        "proxy",
        "locus",
        "allele",
        "frequency",
        "resolution",
        "source_label",
        "source_url",
        "note",
    ]
    out = out.reindex(columns=columns)
    if alleles is not None:
        allele_set = set(alleles)
        out = out[out["allele"].isin(allele_set)].copy()
    return out


def allele_frequency_audit(alleles: Iterable[str]) -> pd.DataFrame:
    """Audit frequency support used for population-coverage calculations.

    The audit keeps regional proxy frequencies separate from published global
    averages. The ``coverage_frequency`` column is the value used by panel
    summaries: population-weighted regional frequency when present, otherwise
    the published global average.
    """
    allele_list = list(dict.fromkeys(alleles))
    rows = allele_frequency_rows(allele_list)
    regional = rows[rows["frequency_scope"] == "regional_proxy"].copy()
    global_df = rows[rows["frequency_scope"] == "published_global"].copy()

    regional_numeric = regional[pd.notna(regional["frequency"])].copy()
    regional_weighted: dict[str, float] = dict.fromkeys(allele_list, 0.0)
    if not regional_numeric.empty:
        per_region = regional_numeric.groupby(["region", "allele"], as_index=False)[
            "frequency"
        ].max()
        total_population = sum(float(population) for population in REGION_POPULATIONS.values())
        for row in per_region.itertuples(index=False):
            population = REGION_POPULATIONS.get(row.region)
            if population is None:
                continue
            regional_weighted[row.allele] += (
                float(row.frequency) * float(population) / total_population
            )

    global_by_allele = global_df[pd.notna(global_df["frequency"])].drop_duplicates(
        subset=["allele"], keep="first"
    )
    global_frequency = {
        str(row.allele): float(row.frequency) for row in global_by_allele.itertuples(index=False)
    }
    global_source_label = {
        str(row.allele): str(row.source_label) for row in global_by_allele.itertuples(index=False)
    }
    global_source_url = {
        str(row.allele): str(row.source_url) for row in global_by_allele.itertuples(index=False)
    }

    audit_rows: list[dict] = []
    for allele in allele_list:
        allele_regional = regional[regional["allele"] == allele]
        regional_frequency = float(regional_weighted.get(allele, 0.0))
        published_global_frequency = global_frequency.get(allele)
        if regional_frequency > 0.0:
            coverage_frequency = regional_frequency
            coverage_source = "regional_weighted"
        elif published_global_frequency is not None:
            coverage_frequency = published_global_frequency
            coverage_source = "published_global"
        else:
            coverage_frequency = 0.0
            coverage_source = "missing"

        audit_rows.append(
            {
                "allele": allele,
                "locus": _allele_locus(allele),
                "regional_weighted_frequency": regional_frequency,
                "published_global_frequency": published_global_frequency,
                "coverage_frequency": coverage_frequency,
                "coverage_frequency_source": coverage_source,
                "numeric_region_count": int(allele_regional["frequency"].notna().sum()),
                "exact_region_count": int((allele_regional["resolution"] == "exact").sum()),
                "proxy_group_region_count": int(
                    (allele_regional["resolution"] == "proxy_group").sum()
                ),
                "qualitative_region_count": int(
                    (allele_regional["resolution"] == "qualitative").sum()
                ),
                "regional_source_labels": "; ".join(
                    sorted(
                        {
                            str(source)
                            for source in allele_regional["source_label"].dropna().tolist()
                        }
                    )
                ),
                "global_source_label": global_source_label.get(allele),
                "global_source_url": global_source_url.get(allele),
            }
        )
    return pd.DataFrame(audit_rows)


def region_names() -> list[str]:
    """Return sorted unique region names from the priority table."""
    return sorted({row["region"] for row in REGION_PRIORITY_ROWS})
