module DeltaQ.Numeric.PDF
  ( module DeltaQ.Numeric.CDF.Types
  , fromPointwisePDF
  )
where

import           Data.List
import qualified Data.Map.Strict as M
import           Data.Ord

import           DeltaQ.Numeric.CDF.Types

-- | Create an empirical CDF from a, potentially improper, point-wise PDF.
fromPointwisePDF :: [(Double, Rational)] -> EmpiricalCDF
fromPointwisePDF [] = makeEmpiricalCDF []
fromPointwisePDF ppdf' =
  checkPdf `seq` makeEmpiricalCDF' Nothing Nothing Nothing cdf
  where
    cdf = M.fromAscList $ scanl1 f' ppdf

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

    f' (_,b) (c,d)
      | d < 0 || d > 1
      = error "fromPointwisePDF: impossible probability mass at point"
      | d' > 1
      = error "fromPointwisePDF: impossible probability mass accumulated"
      | otherwise
      = (c, d')
      where
        d' = b + d
