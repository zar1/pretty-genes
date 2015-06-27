from collections import namedtuple, deque
import matplotlib.pyplot as plt 
import shapely.geometry as sg
import descartes

class Strand:
    """Enum of possibilities for strand"""
    POS, NEG = range(2)

class Marker:
    """Enum of shapes to put at the end of genes

    NONE: no shape, just prints a rectangle, like |=|
    ISOSCELES: gene ends with an isosceles triangle, with a point at the 
        midpoint of the figure. Like |=|>
    RIGHT: gene ends w/ a right traingle w/ a point at the bottom of the 
        figure, like |=|\\
    """

    NONE, ISOSCELES, RIGHT = range(3)
    
    

# Thanks to 
# http://stackoverflow.com/questions/1606436/adding-docstrings-to-namedtuples-in-python 
Gene_ = namedtuple('Gene', ('name', 'strand', 'start', 'end'))
class Gene(Gene_):
    """Class representing a single gene

    Attrs
    -----
    name: str
        name of the gene
    strand: attr of Strand class
    start: int
        Gene's start point in the DNA strand
    end: int
        Gene's end point in the DNA strand
    """    
    pass


class GeneDrawer:

    FormatArgs = namedtuple('FormatArgs', ('marker', 'kwargs'))

    def __init__(self, number_of_lines, **kwargs):
        self.__default_format_kwargs = kwargs
        self.__number_of_lines = number_of_lines
        self.__formats = {}

    def add_format(self, gene_name, marker=Marker.NONE, **kwargs):
        """
        Parameters
        ----------
        gene: str
            Name of the gene to format
        marker: attr of Marker class
            The type of marker for this gene
        kwargs: dict
            kwargs corresponding to kwargs in matplotlib.patches.Rectangle
            (http://matplotlib.org/api/patches_api.html#matplotlib.patches.Rectangle)

        """
        self.__formats[gene_name] = self.FormatArgs(marker, kwargs)

    def __draw_gene(self, gene, torn_left, torn_right, fig):
        gene_name = gene.name
        gene_strand = gene.strand
        gene_start = gene.start
        gene_end = gene.end
        try: 
            format_args = self.__formats[gene_name]
        except KeyError:
            format_args = self.FormatArgs(
                Marker.NONE, 
                self.__default_format_kwargs)
        format_kwargs = format_args.kwargs
        marker = format_args.marker
        upper_left = (gene_start, 0)
        upper_right = (gene_end, 0)
        lower_right = (gene_end, -2)
        lower_left = (gene_start, -2)
        rect = sg.Polygon((upper_left, upper_right, lower_right, lower_left))
        if (marker != Marker.NONE and 
            (not (gene_strand == Strand.NEG and torn_left)) and
            (not (gene_strand == Strand.POS and torn_right))):
            # TODO maybe there's a more elegant way to get our marker width
            x_axis_left, x_axis_right = fig.axes.get_xlim()
            marker_width = (x_axis_right - x_axis_left) * 0.02
            if marker == Marker.RIGHT:
                if gene_strand == Strand.POS:
                    pts = (upper_right, lower_right, (gene_end - marker_width,
                                                      0))
                else:
                    pts = (upper_left, lower_left, (gene_start + marker_width,
                                                    0))
            elif marker == Marker.ISOSCELES:
                if gene_strand == Strand.POS:
                    pts = ((gene_end - marker_width, 0),
                           upper_right,
                           lower_right,
                           (gene_end - marker_width, -2),
                           (gene_end-1, -1))
                else:
                    pts = (upper_left,
                           (gene_start + marker_width, 0),
                           (gene_start+1, -1),
                           (gene_start + marker_width, -2),
                           lower_left)
            marker_cut = sg.Polygon(pts)
            rect = rect.difference(marker_cut)
        # Rather than starting w/ a small rectangle and painting other object on,
        # we start w/ a large rectangle and then cut other objects off
        # Fancy shape manipulation:
        # http://stackoverflow.com/questions/27947054/coloring-intersection-of-circles-patches-in-matplotlib
        #TODO direction, marker, tear, etc.
        #TODO y-values
#        rect = plt.Rectangle((gene_start, -2), gene_end - gene_start,
#                             2, **format_kwargs) 
        fig.add_patch(descartes.PolygonPatch(rect, **format_kwargs))
        #fig.add_patch(descartes.PolygonPatch(marker_cut, fc='g'))

    def __plot_one_line(self, start, end, stream, fig):
            # Make a subplot
            # define our axis
            while stream:
                torn_left, gene = stream.popleft()
                gene_name = gene.name
                gene_strand = gene.strand
                gene_start = gene.start
                gene_end = gene.end
                torn_right = False
                if gene_start >= end:
                    stream.appendleft((torn_left, gene))
                    return
                if gene_end > end:
                    #TODO gene extends accross multiple lines
                    my_part = Gene(gene_name, gene_strand, gene_start, end)
                    their_part = Gene(gene_name, gene_strand, end + 1, 
                                      gene_end)
                    stream.appendleft((True, their_part))
                    gene = my_part
                self.__draw_gene(
                    gene, 
                    torn_left, 
                    torn_right,  
                    fig)
         
    def run(
            self, 
            start, 
            end, 
            genes, 
            figure_kwargs={'figsize': (40, 30), 'dpi': 150}):
        """
        Parameters
        ----------
        start: int
            first gene to display
        end: int
            last gene to display
        genes: list of Gene
            genes to display
        figure_kwargs: dict
            arguments to pass to plt.figure()

        Returns
        -------
        Matplotlib figure     

        """
        fig = plt.figure(**figure_kwargs)
        # the first element of the tuple is whether or not the gene
        # has been torn on the left
        genes_remaining = deque([(False, gene) for gene in genes])
        pairs_per_line = (end - start) / (self.__number_of_lines - 1)
        #TODO what if the right divides the left?
        pairs_last_line = ((end - start) - (pairs_per_line * 
                           self.__number_of_lines))
        this_start = start
        for line in xrange(self.__number_of_lines):
            this_end = this_start + pairs_per_line
            subfig = fig.add_subplot(self.__number_of_lines, 1, line + 1)
            subfig.set_xlim((this_start, this_end))
            subfig.set_ylim((-3, 1))
            subfig.axes.yaxis.set_visible(False)
            subfig.axes.xaxis.tick_top()
            subfig.axes.xaxis.set_label_position('top')
            self.__plot_one_line(
                this_start,
                this_end,
                genes_remaining,
                subfig)
            this_start = this_end + 1
        # TODO treat last line as a special case
#        if pairs_last_line > 0:
#            #TODO cut this axis off earlier
#            subfig = fig.add_subplot(self.__number_of_lines, 1, line + 1)
#            subfig.set_xlim((this_start, this_start + pairs_per_line))
#            subfig.axes.yaxis.set_visible(False)
#            self.__plot_one_line(
#                this_start,
#                this_start + pairs_last_line,
#                genes_remaining,
#                subfig)
        fig.tight_layout()
        return fig
    
if __name__ == '__main__':
    start = 1
    end = 4000000
    genes = [Gene(name, strand, gene_start, gene_end) for 
             name, gene_start, gene_end, strand in
             [('JOHN', 276, 4276, Strand.NEG), ('JOHN', 100000, 104000, Strand.POS), ('ZACH', 1000000, 1004000, Strand.NEG), ('ZACH', 2000000, 2004000, Strand.POS)]]
    d = GeneDrawer(30)
    d.add_format('JOHN', marker=Marker.ISOSCELES, fc='g')
    d.add_format('ZACH', marker=Marker.RIGHT, fc='b')
    plot = d.run(start, end, genes)
    plot.savefig('pretty_genes.eps', format='eps')
