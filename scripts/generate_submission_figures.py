#!/usr/bin/env python3
"""
Generate submission-ready figures for FAdV-4 manuscript
Converts PDF figures to TIFF format at 300 DPI with subfigure labels
"""

import os
import sys
from PIL import Image, ImageDraw, ImageFont
from pdf2image import convert_from_path

# Configuration
DPI = 300
BASE_DIR = "/Users/tgw/Desktop/FADV_new/FAdV4_LimmaAnalysis"
OUTPUT_DIR = os.path.join(BASE_DIR, "submission_files")
FIGURES_DIR = os.path.join(BASE_DIR, "results", "figures")

# Ensure output directory exists
os.makedirs(OUTPUT_DIR, exist_ok=True)

def pdf_to_image(pdf_path, dpi=DPI):
    """Convert PDF to PIL Image at specified DPI"""
    images = convert_from_path(pdf_path, dpi=dpi)
    if images:
        return images[0]
    return None

def add_subfigure_label(img, label):
    """Add subfigure label (A), (B), etc. to top-left corner"""
    draw = ImageDraw.Draw(img)

    # Try to use a bold font, fallback to default
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 48)
    except:
        font = ImageFont.load_default()

    # Draw text with black color
    draw.text((30, 30), label, fill="black", font=font)

    return img

def save_as_tiff(img, filepath):
    """Save image as TIFF with LZW compression"""
    img.save(filepath, format="TIFF", compression="tiff_lzw")
    return filepath

# Define figure configurations
figures_config = {
    "Figure1": {
        "layout": "custom",
        "images": [
            ("placeholder", "A"),  # Will be just label A
            (os.path.join(FIGURES_DIR, "qc/1_PCA_Comparison.pdf"), "B"),
            (os.path.join(FIGURES_DIR, "qc/Statistical_Power_Analysis.pdf"), "C"),
            (os.path.join(FIGURES_DIR, "qc/2_Correlation_Heatmap_Raw.pdf"), "D"),
            (os.path.join(FIGURES_DIR, "qc/3_Correlation_Heatmap_BatchAdjusted.pdf"), "E"),
        ],
        "grid": (2, 3),
        "placeholders": [("A", 0, 0)],  # (label, row, col)
        "descriptions": [
            ("B", 0, 1), ("C", 0, 2),
            ("D", 1, 0), ("E", 1, 1),
        ]
    },
    "Figure2": {
        "layout": "grid",
        "images": [
            (os.path.join(FIGURES_DIR, "comparison/MA_Plot_Combined.pdf"), "A"),
            (os.path.join(FIGURES_DIR, "comparison/Volcano_Plot_Combined.pdf"), "B"),
            (os.path.join(FIGURES_DIR, "Figure_VennDiagram.pdf"), "C"),
        ],
        "grid": (1, 3)
    },
    "Figure3": {
        "layout": "custom_3",
        "images": [
            (os.path.join(FIGURES_DIR, "Figure_Boxplot_Top20_DEGs.pdf"), "A"),
            (os.path.join(FIGURES_DIR, "Enrichment_GO_BP_Bubble.pdf"), "B"),
            (os.path.join(FIGURES_DIR, "Enrichment_KEGG_Bubble.pdf"), "C"),
        ],
    },
    "Figure4": {
        "layout": "custom_4",
        "images": [
            ("placeholder", "A"),
            (os.path.join(FIGURES_DIR, "Figure4C_HubGene_Heatmap.pdf"), "B"),
        ],
    },
    "Figure5": {
        "layout": "single",
        "images": [
            (os.path.join(FIGURES_DIR, "WGCNA_Module_Trait_Heatmap.pdf"), ""),
        ],
    },
    "Supplementary_Figure_S1": {
        "layout": "single",
        "images": [
            (os.path.join(FIGURES_DIR, "WGCNA_SoftThreshold.pdf"), ""),
        ],
    },
    "Supplementary_Figure_S2": {
        "layout": "single",
        "images": [
            (os.path.join(FIGURES_DIR, "WGCNA_SampleClustering.pdf"), ""),
        ],
    },
    "Supplementary_Figure_S3": {
        "layout": "vertical_2",
        "images": [
            (os.path.join(FIGURES_DIR, "WGCNA_Module_Dendrogram.pdf"), "A"),
            (os.path.join(FIGURES_DIR, "WGCNA_Eigengene_Adjacency.pdf"), "B"),
        ],
    },
}

def create_figure1():
    """Create Figure 1: 2x3 grid with placeholder A"""
    print("Creating Figure1...")

    # Convert all PDFs to images
    images = {}
    labels = {}

    # B: PCA
    img_b = pdf_to_image(os.path.join(FIGURES_DIR, "qc/1_PCA_Comparison.pdf"))
    if img_b:
        images["B"] = img_b
        print(f"  - B (PCA): {img_b.size}")

    # C: Power
    img_c = pdf_to_image(os.path.join(FIGURES_DIR, "qc/Statistical_Power_Analysis.pdf"))
    if img_c:
        images["C"] = img_c
        print(f"  - C (Power): {img_c.size}")

    # D: Correlation Raw
    img_d = pdf_to_image(os.path.join(FIGURES_DIR, "qc/2_Correlation_Heatmap_Raw.pdf"))
    if img_d:
        images["D"] = img_d
        print(f"  - D (Raw): {img_d.size}")

    # E: Correlation Adjusted
    img_e = pdf_to_image(os.path.join(FIGURES_DIR, "qc/3_Correlation_Heatmap_BatchAdjusted.pdf"))
    if img_e:
        images["E"] = img_e
        print(f"  - E (Adjusted): {img_e.size}")

    # Create grid: 2 rows, 3 columns
    # Layout:
    # [ A (placeholder) ] [ B PCA ] [ C Power ]
    # [ D Raw ]           [ E Adj ] [ empty   ]

    # Use B as reference size
    ref_width = images["B"].width
    ref_height = images["B"].height

    # Calculate cell size (add spacing)
    spacing = 30
    cell_width = ref_width
    cell_height = ref_height

    # Create canvas: 3 cells wide, 2 cells tall
    canvas_width = cell_width * 3 + spacing * 4
    canvas_height = cell_height * 2 + spacing * 3
    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")
    draw = ImageDraw.Draw(canvas)

    # Try font
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
    except:
        font = ImageFont.load_default()

    # Place images
    # Row 0: A(empty), B, C
    # Row 1: D, E, empty

    # A: Placeholder - just draw "A" label
    a_x = spacing
    a_y = spacing
    draw.text((a_x + 20, a_y + 20), "A", fill="black", font=font)

    # B: PCA at (col 1, row 0)
    b_x = cell_width + spacing * 2
    b_y = spacing
    if "B" in images:
        canvas.paste(images["B"], (b_x, b_y))
        draw = ImageDraw.Draw(canvas)
        draw.text((b_x + 30, b_y + 30), "B", fill="black", font=font)

    # C: Power at (col 2, row 0)
    c_x = cell_width * 2 + spacing * 3
    c_y = spacing
    if "C" in images:
        canvas.paste(images["C"], (c_x, c_y))
        draw = ImageDraw.Draw(canvas)
        draw.text((c_x + 30, c_y + 30), "C", fill="black", font=font)

    # D: Raw at (col 0, row 1)
    d_x = spacing
    d_y = cell_height + spacing * 2
    if "D" in images:
        canvas.paste(images["D"], (d_x, d_y))
        draw = ImageDraw.Draw(canvas)
        draw.text((d_x + 30, d_y + 30), "D", fill="black", font=font)

    # E: Adjusted at (col 1, row 1)
    e_x = cell_width + spacing * 2
    e_y = cell_height + spacing * 2
    if "E" in images:
        canvas.paste(images["E"], (e_x, e_y))
        draw = ImageDraw.Draw(canvas)
        draw.text((e_x + 30, e_y + 30), "E", fill="black", font=font)

    output_path = os.path.join(OUTPUT_DIR, "Figure1.tiff")
    save_as_tiff(canvas, output_path)
    print(f"  Saved: {output_path} ({canvas.size})")

    return canvas

def create_figure2():
    """Create Figure 2: 1x3 grid"""
    print("Creating Figure2...")

    pdfs = [
        (os.path.join(FIGURES_DIR, "comparison/MA_Plot_Combined.pdf"), "A"),
        (os.path.join(FIGURES_DIR, "comparison/Volcano_Plot_Combined.pdf"), "B"),
        (os.path.join(FIGURES_DIR, "Figure_VennDiagram.pdf"), "C"),
    ]

    images = []
    for pdf_path, label in pdfs:
        img = pdf_to_image(pdf_path)
        if img:
            print(f"  - {label}: {img.size}")
            images.append((img, label))

    if not images:
        print("  ERROR: No images converted!")
        return

    # Get reference size from first image
    ref_width = images[0][0].width
    ref_height = images[0][0].height

    # Resize all to same height
    target_height = ref_height
    resized = []
    for img, label in images:
        if img.height != target_height:
            ratio = target_height / img.height
            new_width = int(img.width * ratio)
            img = img.resize((new_width, target_height), Image.LANCZOS)
        resized.append((img, label))

    # Create canvas: 3 columns, 1 row
    spacing = 30
    cell_width = max(img.width for img, _ in resized)
    canvas_width = cell_width * 3 + spacing * 4
    canvas_height = target_height + spacing * 2

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
    except:
        font = ImageFont.load_default()

    for i, (img, label) in enumerate(resized):
        x = spacing + i * (cell_width + spacing)
        y = spacing
        canvas.paste(img, (x, y))
        draw = ImageDraw.Draw(canvas)
        draw.text((x + 30, y + 30), label, fill="black", font=font)

    output_path = os.path.join(OUTPUT_DIR, "Figure2.tiff")
    save_as_tiff(canvas, output_path)
    print(f"  Saved: {output_path} ({canvas.size})")

    return canvas

def create_figure3():
    """Create Figure 3: A on top (full width), B and C below side by side"""
    print("Creating Figure3...")

    # A: Top 20 DEGs boxplot
    img_a = pdf_to_image(os.path.join(FIGURES_DIR, "Figure_Boxplot_Top20_DEGs.pdf"))
    # B: GO BP bubble
    img_b = pdf_to_image(os.path.join(FIGURES_DIR, "Enrichment_GO_BP_Bubble.pdf"))
    # C: KEGG bubble
    img_c = pdf_to_image(os.path.join(FIGURES_DIR, "Enrichment_KEGG_Bubble.pdf"))

    if not all([img_a, img_b, img_c]):
        print("  ERROR: Missing images!")
        return

    print(f"  - A: {img_a.size}")
    print(f"  - B: {img_b.size}")
    print(f"  - C: {img_c.size}")

    spacing = 30

    # Top row: A takes full width
    # Bottom row: B and C side by side

    # Use A as reference width
    ref_width = img_a.width

    # Resize B and C to match A's width/2
    half_width = (ref_width - spacing) // 2

    # Resize images proportionally
    def resize_keep_aspect(img, target_width):
        ratio = target_width / img.width
        new_height = int(img.height * ratio)
        return img.resize((target_width, new_height), Image.LANCZOS)

    # For A: keep original width, scale height proportionally
    target_height_a = int(img_a.height * (ref_width / img_a.width))
    img_a_resized = img_a.resize((ref_width, target_height_a), Image.LANCZOS)

    # For B and C: half width
    img_b_resized = resize_keep_aspect(img_b, half_width)
    img_c_resized = resize_keep_aspect(img_c, half_width)

    # Calculate canvas size
    canvas_width = ref_width + spacing * 2
    canvas_height = img_a_resized.height + max(img_b_resized.height, img_c_resized.height) + spacing * 3

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
    except:
        font = ImageFont.load_default()

    # Place A at top
    x_a = spacing
    y_a = spacing
    canvas.paste(img_a_resized, (x_a, y_a))
    draw = ImageDraw.Draw(canvas)
    draw.text((x_a + 30, y_a + 30), "A", fill="black", font=font)

    # Place B and C at bottom
    y_bottom = img_a_resized.height + spacing * 2

    x_b = spacing
    canvas.paste(img_b_resized, (x_b, y_bottom))
    draw = ImageDraw.Draw(canvas)
    draw.text((x_b + 30, y_bottom + 30), "B", fill="black", font=font)

    x_c = spacing + half_width + spacing
    canvas.paste(img_c_resized, (x_c, y_bottom))
    draw = ImageDraw.Draw(canvas)
    draw.text((x_c + 30, y_bottom + 30), "C", fill="black", font=font)

    output_path = os.path.join(OUTPUT_DIR, "Figure3.tiff")
    save_as_tiff(canvas, output_path)
    print(f"  Saved: {output_path} ({canvas.size})")

    return canvas

def create_figure4():
    """Create Figure 4: Placeholder A and Hub gene heatmap B"""
    print("Creating Figure4...")

    # B: Hub gene heatmap
    img_b = pdf_to_image(os.path.join(FIGURES_DIR, "Figure4C_HubGene_Heatmap.pdf"))

    if not img_b:
        print("  ERROR: Hub gene heatmap not found!")
        return

    print(f"  - B: {img_b.size}")

    spacing = 30
    cell_width = img_b.width
    cell_height = img_b.height

    # Create canvas: 2 columns, 1 row
    canvas_width = cell_width * 2 + spacing * 3
    canvas_height = cell_height + spacing * 2

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
    except:
        font = ImageFont.load_default()

    # A: Placeholder at column 0
    draw = ImageDraw.Draw(canvas)
    draw.text((spacing + 20, spacing + 20), "A", fill="black", font=font)

    # B: Hub gene heatmap at column 1
    x_b = cell_width + spacing * 2
    y_b = spacing
    canvas.paste(img_b, (x_b, y_b))
    draw = ImageDraw.Draw(canvas)
    draw.text((x_b + 30, y_b + 30), "B", fill="black", font=font)

    output_path = os.path.join(OUTPUT_DIR, "Figure4.tiff")
    save_as_tiff(canvas, output_path)
    print(f"  Saved: {output_path} ({canvas.size})")

    return canvas

def create_single_figure(name, pdf_path):
    """Create a single figure from PDF"""
    print(f"Creating {name}...")

    img = pdf_to_image(pdf_path)
    if not img:
        print(f"  ERROR: {pdf_path} not found!")
        return

    print(f"  - Source: {img.size}")

    # Add label if needed
    draw = ImageDraw.Draw(img)
    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 48)
    except:
        font = ImageFont.load_default()
    draw.text((30, 30), "", fill="black", font=font)

    output_path = os.path.join(OUTPUT_DIR, f"{name}.tiff")
    save_as_tiff(img, output_path)
    print(f"  Saved: {output_path} ({img.size})")

    return img

def create_supplementary_figure_s3():
    """Create Supplementary Figure S3: Two images stacked vertically"""
    print("Creating Supplementary_Figure_S3...")

    img_a = pdf_to_image(os.path.join(FIGURES_DIR, "WGCNA_Module_Dendrogram.pdf"))
    img_b = pdf_to_image(os.path.join(FIGURES_DIR, "WGCNA_Eigengene_Adjacency.pdf"))

    if not all([img_a, img_b]):
        print("  ERROR: Missing images!")
        return

    print(f"  - A (Dendrogram): {img_a.size}")
    print(f"  - B (Adjacency): {img_b.size}")

    spacing = 30

    # Resize to same width
    target_width = max(img_a.width, img_b.width)

    def resize_to_width(img, target_w):
        ratio = target_w / img.width
        new_height = int(img.height * ratio)
        return img.resize((target_w, new_height), Image.LANCZOS)

    img_a_resized = resize_to_width(img_a, target_width)
    img_b_resized = resize_to_width(img_b, target_width)

    # Canvas
    canvas_width = target_width + spacing * 2
    canvas_height = img_a_resized.height + img_b_resized.height + spacing * 3

    canvas = Image.new("RGB", (canvas_width, canvas_height), "white")

    try:
        font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 72)
    except:
        font = ImageFont.load_default()

    # Place A at top
    x = spacing
    y = spacing
    canvas.paste(img_a_resized, (x, y))
    draw = ImageDraw.Draw(canvas)
    draw.text((x + 30, y + 30), "A", fill="black", font=font)

    # Place B at bottom
    y_b = img_a_resized.height + spacing * 2
    canvas.paste(img_b_resized, (x, y_b))
    draw = ImageDraw.Draw(canvas)
    draw.text((x + 30, y_b + 30), "B", fill="black", font=font)

    output_path = os.path.join(OUTPUT_DIR, "Supplementary_Figure_S3.tiff")
    save_as_tiff(canvas, output_path)
    print(f"  Saved: {output_path} ({canvas.size})")

    return canvas

def main():
    print("=" * 60)
    print("Generating submission-ready figures")
    print("=" * 60)
    print(f"Output directory: {OUTPUT_DIR}")
    print(f"DPI: {DPI}")
    print()

    # Create all figures
    create_figure1()
    print()

    create_figure2()
    print()

    create_figure3()
    print()

    create_figure4()
    print()

    # Single figures
    create_single_figure("Figure5", os.path.join(FIGURES_DIR, "WGCNA_Module_Trait_Heatmap.pdf"))
    print()

    create_single_figure("Supplementary_Figure_S1", os.path.join(FIGURES_DIR, "WGCNA_SoftThreshold.pdf"))
    print()

    create_single_figure("Supplementary_Figure_S2", os.path.join(FIGURES_DIR, "WGCNA_SampleClustering.pdf"))
    print()

    create_supplementary_figure_s3()
    print()

    print("=" * 60)
    print("All figures generated successfully!")
    print("=" * 60)

    # List output files
    print("\nOutput files:")
    for f in sorted(os.listdir(OUTPUT_DIR)):
        if f.endswith(".tiff"):
            fpath = os.path.join(OUTPUT_DIR, f)
            size = os.path.getsize(fpath)
            print(f"  {f} ({size / 1024:.1f} KB)")

if __name__ == "__main__":
    main()
