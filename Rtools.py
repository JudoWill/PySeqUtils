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


def convert_columns_to_factors(rpy_dataframe, col_refs):
    """Converts a set of columns in a RPY dataframe into Factors.

      col_refs -- [(colA, refA), (colB, refB), ...]
    """

    relevel_fun = robjects.r("""relevel_obj <- function(inputframe, column, ref)
{
inputframe[,column] <- relevel(factor(inputframe[,column]), ref = ref)
return(inputframe)
}""")

    for col, ref in col_refs:
        rpy_dataframe = relevel_fun(rpy_dataframe, col, ref)

    return rpy_dataframe


def R_linear_mixed_effects_model(rpy_dataframe, model_eqn, random_eqn):
    """Runs the lme R model using the input RPY dataframe, model equation and random effects.
    """

    _ = robjects.r('require(nlme)')
    lme_model_func = robjects.r("""lmefunc <- function(data, eqn, reqn){
lm1 <- lme(eqn, random = reqn, data = data, na.action = na.omit)
return(summary(lm1))
}""")

    summary = lme_model_func(rpy_dataframe, model_eqn, random_eqn)
    outdata = {}
    for name in summary.names:
        try:
            outdata[name] = com.convert_robj(summary.rx2(name))
        except:
            pass

    return outdata
