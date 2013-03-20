__author__ = 'will'

#from pandas import DataFrame
import pandas.rpy.common as com
import rpy2.robjects as robjects


def quantile_norm_with_R(input_df):
    """Uses R normalize.quantiles to normalize a DataFrame.
    """

    robjects.r('require(preprocessCore)')
    R_norm_func = robjects.r("""quantnorm <- function(inputmatrix)
{
y<-normalize.quantiles(inputmatrix)
return(y)
}""")

    R_matrix = com.convert_to_r_matrix(input_df)

    normed_matrix = R_norm_func(R_matrix)
    normed_df = com.convert_robj(normed_matrix)

    normed_df.index = input_df.index
    normed_df.columns = input_df.columns

    return normed_df