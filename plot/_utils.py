from matplotlib.colors import to_hex,LinearSegmentedColormap
import matplotlib.pyplot as plt


from typing import Union, Literal, List, Dict, Optional



def generate_gradient_colors(basecolor, ncolor = 10, light_factor = 0.3, reverse = False):
    if ncolor != None:
        ncolor += ncolor//5
        cmap = LinearSegmentedColormap.from_list('custom_gradient', [basecolor, 'white'], N = ncolor)
        gradient_colors_hex = [to_hex(cmap(i/ncolor)) for i in range(ncolor)]
        return gradient_colors_hex
    else:
        from matplotlib.colors import to_rgba
        base_rgba = to_rgba(basecolor)
        light_color = tuple(light_factor * comp + (1 - light_factor) for comp in base_rgba[:3]) 
        if not reverse:
            cmap = LinearSegmentedColormap.from_list('custom_gradient', [basecolor, light_color])
        else:
            cmap = LinearSegmentedColormap.from_list('custom_gradient', [light_color, basecolor])
        return cmap


def get_colors(
    total_for_plot: Union[int,List], 
    colors: Optional[Union[str, List[str], Dict[str, str]]] = None,
):
    if isinstance(total_for_plot, int):
        total_for_plot = [1] * total_for_plot
        
    if colors == None:
        # colors = plt.get_cmap('tab10')(range(len(total_for_plot)))
        cmap = plt.get_cmap('tab20')
        colors = [cmap(i % cmap.N) for i in range(len(total_for_plot))]
        colors = [to_hex(c) for c in colors]
        
    elif isinstance(colors,str):
        if len(total_for_plot) == 1:
            return [colors]
        try:
            cmap = plt.get_cmap(colors)
            colors = [cmap(i % cmap.N) for i in range(len(total_for_plot))]
            colors = [to_hex(c) for c in colors]
        
        except:
            colors = generate_gradient_colors(basecolor = colors, ncolor = len(total_for_plot))
        
    elif isinstance(colors,list):
        if len(colors) < len(total_for_plot):
            raise ValueError("Lenth not match bewteen colors and chains")
            
    elif isinstance(colors,dict):
        colors = [ colors[plot] for plot in total_for_plot]
    return colors

def adjust_font_size(ax):
    fig = ax.figure
    dpi = fig.dpi
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    base_size = 12
    font_size = max(width, height) * base_size / 4.5
    
    return font_size

def _get_color_list(k):
    color_pal1 = ["#70f3ff",  "#44cef6",  "#3eede7",  "#1685a9",  "#177cb0",  "#065279",  "#003472",
  "#4b5cc4",  "#a1afc9",  "#2e4e7e",  "#3b2e7e",  "#4a4266",  "#426666",  "#425066",  "#574266", 
  "#8d4bbb",  "#815463",  "#815476",  "#4c221b",  "#003371",  "#56004f",  "#801dae",  "#4c8dae",
  "#b0a4e3",  "#cca4e3",  "#edd1d8",  "#e4c6d0",  "#ff461f",  "#ff2d51",  "#f36838",  "#ed5736",
  "#ff4777",  "#f00056",  "#ffb3a7",  "#f47983",  "#db5a6b",  "#c93756",  "#f9906f",  "#f05654",
  "#ff2121",  "#f20c00",  "#8c4356",  "#c83c23",  "#9d2933",  "#ff4c00",  "#ff4e20",  "#f35336",
  "#dc3023",  "#ff3300",  "#cb3a56",  "#a98175",  "#b36d61",  "#ef7a82",  "#ff0097",  "#c32136",
  "#be002f",  "#c91f37",  "#bf242a",  "#c3272b",  "#9d2933",  "#60281e",  "#622a1d",  "#bce672",
  "#c9dd22",  "#bddd22",  "#afdd22",  "#a3d900",  "#9ed900",  "#9ed048",  "#96ce54",  "#00bc12",
  "#0eb83a",  "#0eb83a",  "#0aa344",  "#16a951",  "#21a675",  "#057748",  "#0c8918",  "#00e500",
  "#40de5a",  "#00e079",  "#00e09e",  "#3de1ad",  "#2add9c",  "#2edfa3",  "#7fecad",  "#a4e2c6",
  "#7bcfa6",  "#1bd1a5",  "#48c0a3",  "#549688",  "#789262",  "#758a99",  "#50616d",  "#424c50",
  "#41555d",  "#eaff56",  "#fff143",  "#faff72",  "#ffa631",  "#ffa400",  "#fa8c35",  "#ff8c31",
  "#ff8936",  "#ff7500",  "#ffb61e",  "#ffc773",  "#ffc64b",  "#f2be45",  "#f0c239",  "#e9bb1d", 
  "#d9b611",  "#eacd76",  "#eedeb0",  "#d3b17d",  "#e29c45",  "#a78e44",  "#c89b40",  "#ae7000", 
  "#ca6924",  "#b25d25",  "#b35c44",  "#9b4400",  "#9c5333",  "#a88462",  "#896c39",  "#827100",  
  "#6e511e",  "#7c4b00",  "#955539",  "#845a33"]
    color_pal1 = list(set(color_pal1)) * 100
    def generate_color_list(colors, k, seed=61):
        random.seed(seed)
        color_list = random.sample(colors, k)
        while any(color_list[i] == color_list[i + 1] for i in range(len(color_list) - 1)):
            random.shuffle(color_list)
        return color_list

    import random
    random.seed(61)

    color_lists = generate_color_list(color_pal1, k)

    return color_lists