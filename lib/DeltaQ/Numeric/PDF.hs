module DeltaQ.Numeric.PDF
  ( module DeltaQ.Numeric.CDF.Types
  , fromPointwisePDF
  )
where

import Data.List
import Data.Ord

import DeltaQ.Numeric.CDF.Types

-- | Create an empirical CDF from a, potentially improper, point-wise PDF.
fromPointwisePDF :: [(Double, Rational)] -> EmpiricalCDF
fromPointwisePDF [] = makeEmpiricalCDF []
fromPointwisePDF ppdf' =
  checkPdf `seq` construct cdf
  where
    ppdf = sortBy (comparing fst) ppdf'

    checkPdf
      | b < 0 || b > 1
      = error "fromPointwisePDF: impossible probability mass at point"
      | a < 0
      = error "fromPointwisePDF: attempted causality breach detected"
      | otherwise
        = ()
      where
        (a,b) = head ppdf

    cdf = scanl1 f' ppdf

    f' (_,b) (c,d)
      | d < 0 || d > 1
      = error "fromPointwisePDF: impossible probability mass at point"
      | d' > 1
      = error "fromPointwisePDF: impossible probability mass accumulated"
      | otherwise
      = (c, d')
      where
        d' = b + d

    -- There is definitely room for computational efficiency improvement here.
    -- Multiple passes over the data just make the intentions clearer. NOTE: we
    -- punted the empty list case, so there is (at least one) data point
    -- present.
    construct :: [(Double, Rational)] -> EmpiricalCDF
    construct _ =   result
      where
        result = EmpiricalCDF
          { _ecdf   = f
          , _ecdfInverse = g
          , _ecdfPdf     = undefined
          , _ecdfMass    = undefined
          , _ecdfSamples = undefined
          , _ecdfMin     = undefined
          , _ecdfMax     = undefined
          , _ecdfMean    = undefined
          , _ecdfVariance = undefined
          }
        f x
          | x <= _ecdfMin result = 0
          | x >= _ecdfMax result = _ecdfMass result
          | otherwise = undefined
        g x
          | x <  0    = error "makeEmpiricalCDF: negative probability"
          | x >  1    = error "makeEmpiricalCDF: > unit probability"
          | x >  _ecdfMass result    = Nothing
          | otherwise = undefined
