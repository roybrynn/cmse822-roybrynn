import csv
import matplotlib.pyplot as plt

# Use the strong scaling CSV file
csv_path = "../agoge/scaling_logs/zone_updates_strong_nonblocking_periodic.csv"  # adjust path if needed

# Read CSV using the built-in csv module and group data by baseline (zones per process).
data = {}
with open(csv_path, newline='') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        # Use "Baseline" as the grouping key; convert to int for proper sorting.
        baseline = int(row["Baseline"].strip())
        procs = int(row["Processors"].strip())
        zone_updates = float(row["ZoneUpdates(M)"].strip())
        if baseline not in data:
            data[baseline] = []
        data[baseline].append((procs, zone_updates))

# Sort data points for each baseline group by processor count.
for baseline in data:
    data[baseline].sort(key=lambda tup: tup[0])

# Set up a marker cycle for different baselines
markers = ['o', 's', '^', 'D', 'v', 'P', '*', 'X']

# ----------------- First Plot: Raw Zone Updates -----------------
plt.figure(figsize=(8, 6))
unique_baselines = sorted(data.keys())
for i, baseline in enumerate(unique_baselines):
    procs, updates = zip(*data[baseline])
    linestyle = '-' + markers[i % len(markers)]
    plt.plot(procs, updates, linestyle, label=f"Baseline={baseline}", markersize=8)
plt.xlabel("Processor Count", fontsize=12)
plt.ylabel("Zone Updates per Second (M)", fontsize=12)
plt.title("Strong Scaling: Zone Updates vs. Processor Count", fontsize=14)
plt.legend(title="Zones per Process", fontsize=10)
plt.loglog()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("strong_performance_plot.png", dpi=300)
# plt.show()

# ----------------- Second Plot: Normalized by Processor Count -----------------
plt.figure(figsize=(8, 6))
normalized_data = {}  # Save normalized values for efficiency plot later.
for i, baseline in enumerate(unique_baselines):
    procs, updates = zip(*data[baseline])
    norm_updates = [u / p for p, u in zip(procs, updates)]
    normalized_data[baseline] = (procs, norm_updates)
    linestyle = '-' + markers[i % len(markers)]
    plt.plot(procs, norm_updates, linestyle, label=f"Baseline={baseline}", markersize=8)
plt.xlabel("Processor Count", fontsize=12)
plt.ylabel("Zone Updates per Core-second (M)", fontsize=12)
plt.title("strong Scaling: Normalized Zone Updates per Core-second", fontsize=14)
plt.legend(title="Zones per Process", fontsize=10)
plt.loglog()
plt.grid(True, linestyle="--", alpha=0.6)
plt.tight_layout()
plt.savefig("strong_normalized_performance_plot.png", dpi=300)
# plt.show()

# ----------------- Third Plot: Parallel Efficiency -----------------
plt.figure(figsize=(8, 6))
for i, baseline in enumerate(unique_baselines):
    procs, norm_updates = normalized_data[baseline]
    # Use the normalized value at the lowest processor count as baseline efficiency.
    base = norm_updates[0]
    efficiency = [val / base for val in norm_updates]
    linestyle = '-' + markers[i % len(markers)]
    plt.plot(procs, efficiency, linestyle, label=f"Baseline={baseline}", markersize=8)
plt.xlabel("Processor Count", fontsize=12)
plt.ylabel("Parallel Efficiency", fontsize=12)
plt.title("strong Scaling: Parallel Efficiency vs. Processor Count", fontsize=14)
plt.legend(title="Zones per Process", fontsize=10)
plt.grid(True, linestyle="--", alpha=0.6)
plt.semilogx()  # x-axis log scale; linear y-axis
plt.tight_layout()
plt.savefig("strong_parallel_efficiency_plot.png", dpi=300)
plt.show()
