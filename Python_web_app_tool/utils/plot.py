import matplotlib.pyplot as plt
import io
import base64
from markupsafe import Markup

def plot_distribution(counts, title, x_label='Nucleotide', y_label='Percentage'):
    """Create a bar chart and return it as an HTML-safe image."""
    plt.figure(figsize=(5, 5))
    plt.bar(counts.keys(), counts.values())
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(title)
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight')
    buf.seek(0)
    image_base64 = base64.b64encode(buf.read()).decode('utf-8')
    buf.close()
    plt.close()
    return Markup(f'<img src="data:image/png;base64,{image_base64}" />')

plot_nucleotide_distribution = lambda counts: plot_distribution(counts, "Nucleotide Composition")
plot_amino_acid_frequencies = lambda counts: plot_distribution(
    dict(sorted(counts.items(), key=lambda x: x[1], reverse=True)[:20]),
    "Top 20 Amino Acid Frequencies", x_label='Amino Acids', y_label='Frequency')