from collections import namedtuple, deque
import matplotlib.pyplot as plt 

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
        #TODO direction, marker, tear, etc.
        #TODO y-values
        rect = plt.Rectangle((gene_start, -2), gene_end - gene_start,
                             2, **format_kwargs) 
        fig.add_patch(rect)

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
                if gene_start > end:
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
         
    def run(self, start, end, genes):
        """
        Parameters
        ----------
        start: int
            first gene to display
        end: int
            last gene to display
        genes: list of Gene
            genes to display

        Returns
        -------
        Matplotlib figure     

        """
        fig = plt.figure()
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
            subfig.set_ylim((-3, 0))
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
        return fig
    
if __name__ == '__main__':
    start = 1
    end = 3000
    genes = [Gene('JOHNGENE', Strand.NEG, gene_start, gene_end) for 
             gene_start, gene_end in
             [(276, 518), (732, 1000), (1500, 2000), (2200, 3000)]]
    d = GeneDrawer(5)
    plot = d.run(start, end, genes)
    plot.show()
    import pdb; pdb.set_trace()
