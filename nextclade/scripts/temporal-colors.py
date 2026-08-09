"""
Custom script to read dates (YYYY-MM-DD format) from a TSV file,
sort them, and create a custom ordering & colour map for use
within an auspice-config JSON
"""

import argparse
import csv
import json
import math
from collections import Counter

rainbow = ["#511EA8", "#4928B4", "#4334BF", "#4041C7", "#3F50CC", "#3F5ED0", "#416CCE", "#4379CD", "#4784C7", "#4B8FC1", "#5098B9", "#56A0AF", "#5CA7A4", "#63AC99", "#6BB18E", "#73B583", "#7CB878", "#86BB6E", "#90BC65", "#9ABD5C", "#A4BE56", "#AFBD4F", "#B9BC4A", "#C2BA46", "#CCB742", "#D3B240", "#DAAC3D", "#DFA43B", "#E39B39", "#E68F36", "#E68234", "#E67431", "#E4632E", "#E1512A", "#DF4027", "#DC2F24"]
missing_color = "#ADB1B3" # same as Auspice <https://github.com/nextstrain/auspice/blob/e9e910bb9000e17173aac58d7399893003cc8af1/src/util/colorScale.ts#L17>
missing_date = "XXXX-XX-XX"

def assign_colors(dates, k):
    """
    Project the observed (numeric) date range into the rainbow colour range using an
    exponential projection so that newer dates are more spread out.
    k controls curvature: k=0 is linear, higher k compresses older dates more."""
    n = len(rainbow)
    min_val = dates[0]['numeric']
    max_val = dates[-1]['numeric']
    for d in dates:
        t = (d['numeric'] - min_val) / (max_val - min_val)
        if k == 0:
            d['idx'] = math.floor(t * (n - 1))
        else:
            d['idx'] = math.floor((math.exp(k * t) - 1) / (math.expm1(k)) * (n - 1))
        d['color'] = rainbow[d['idx']]
        

def numeric(d: str) -> float:
    """convert YYYY-MM-DD input to numeric output. Ignore leap-year complexity."""
    year, month, day = d.split('-')
    y = int(year)
    if month == 'XX':
        return y + 0.5
    elif day == 'XX':
        m = int(month)
        days_in_month = [31,28,31,30,31,30,31,31,30,31,30,31][m-1]
        start = sum([31,28,31,30,31,30,31,31,30,31,30,31][:m-1])
        return y + (start + days_in_month / 2) / 365
    else:
        m = int(month)
        dd = int(day)
        day_of_year = sum([31,28,31,30,31,30,31,31,30,31,30,31][:m-1]) + dd
        return y + day_of_year / 365

def auspice_json(dates, key):
    # Note: if there's LOTS of dates we can control the legend entries here. See https://github.com/nextstrain/augur/blob/5b96f9ea89a602711fa2fd7d3068762b926af5db/augur/data/schema-auspice-config-v2.json#L54
    coloring = {
        "key": key,
        "title": "Collection Date",
        "type": "categorical",
        "scale": [ [d["value"], d["color"] ] for d in dates]
    }
    return coloring
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--metadata', required=True, type=str, help="Input metadata TSV")
    parser.add_argument('--key', default="date", type=str, help="TSV column to use")
    parser.add_argument('--output', required=True, type=str, help="JSON output path")
    args = parser.parse_args()

    with open(args.metadata, newline='') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        values = [row[args.key] for row in reader if row.get(args.key)]

    counts = Counter(values)

    dates = sorted(
        [{"numeric": numeric(d), "count": counts[d], "value": d}
         for d in set(values) if d != missing_date],
        key=lambda x: x["numeric"]
    )

    assign_colors(dates, 2.0)

    if missing_date in values:
        dates.insert(0, {"count": counts[missing_date], "value": missing_date, "color": missing_color})

    # for d in dates:
    #     print(d)

    coloring = auspice_json(dates, args.key)
    with open(args.output, 'w') as fh:
        json.dump(coloring, fh, indent=2)
