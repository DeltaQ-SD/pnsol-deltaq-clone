{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE TemplateHaskell #-}
{-# LANGUAGE LambdaCase #-}
module DeltaQ.Numeric.CDF.Types
where

import           Control.Lens
import           Control.Monad.State.Strict
import           Data.Default
import qualified Data.Map.Strict as M
import           Data.Maybe
import           Numeric.IEEE
import           Text.Show.Functions ()

-- | An empirical CDF for improper random variable.
data EmpiricalCDF = EmpiricalCDF
  { _ecdf         :: !(Double -> Rational)
    -- ^ The empirical CDF, may not reach 1. It is undefined if the
    --   `ecdfMass` is zero.
  , _ecdfInverse  :: !(Rational -> Maybe Double)
    -- ^ Returns values over [0,`ecdfMass`] and `Nothing` between
    --   (`ecdMass`, 1]. The inverse of the above.
  , _ecdfPdf       :: (Double -> ((Double, Double), Rational))
    -- ^ The empirical PDF derived from the CDF. Note that this is non-strictly
    --   evaluated (to defer cost if not used). The resulting function returns a
    --   probability mass density over the half-open interval [a,b) that
    --   contains the given point.
  , _ecdfMass     :: !Rational
    -- ^ The tangible propabilty mass.
  , _ecdfSamples  :: !Int
    -- ^ The total number of samples taken for the construction of
    --   `EmpiricalCDF`.
  , _ecdfMin      :: !Double
    -- ^ The minimum sampled value. Value is +∞ if `ecdfMass` is zero.
  , _ecdfMax      :: !Double
    -- ^ The minimum sampled value. Value is -∞ if `ecdfMass` is zero.
  , _ecdfMean     :: !Double
    -- ^ The mean of the sampled tangible distribution. Zero if
    --   `ecdfMass` is zero.
  , _ecdfVariance :: !Double
    -- ^ The variance of the sampled tangible distribution. Is NaN if
    --   the number samples is two or less.
  }
    deriving Show
makeLenses ''EmpiricalCDF

-- | Given a list of IRV samples construct an empirical improper CDF
--   along with some additional statistics.
makeEmpiricalCDF :: [Maybe Double] -> EmpiricalCDF
makeEmpiricalCDF ys = run (step ys >> finalise)
  where
    -- execute the evalation loop
    run  = (flip evalState) def
    -- loop over the input stream
    step [] = return ()
    step (Nothing:zs)
      = _2 += 1 >> step zs -- just increment the total count
    step ((Just z):zs) = do
      -- insert in the value to occurance map
      _1 %= M.insertWith (+) z 1
      -- increment the total count
      _2 += 1
      -- increment the tangible mass count
      _3 += (1 :: Int)
      -- step the online mean and variance algorithm
      do n      <- fmap fromIntegral $ use _3
         delta  <- fmap (z -) $ use _4
         _4     += delta / n
         delta2 <- fmap (z -) $ use _4
         _5     += delta * delta2
      -- and loop
      step zs
    finalise = do
      (!m',!n,!t,!a,!b) <- get
          -- finalise the variance
      let v | t > 2     = b / fromIntegral (t - 1)
            | otherwise = nan
          -- normalise the occurance map into cumulative probability
          m = normalise n m'
      return $ makeEmpiricalCDF' (Just n) (Just a) (Just v) m
    normalise :: Int -> M.Map Double Int -> M.Map Double Rational
    normalise n
      = snd
      . M.mapAccum (\a b -> let c = a + toRational b in (c, c / toRational n)) 0

-- | Generate an EmpiricalCDF from a normalized CDF map. As initial construction
--   can more easily do some calculations, permit those to be injected.
makeEmpiricalCDF' :: Maybe Int
                  -- ^ The number of samples (Nothing implies size of supplied map)
                  -> Maybe Double
                  -- ^ The calculated arithmetic mean (`Nothing` implies it is
                  -- calculated from the supplied map)
                  -> Maybe Double
                  -- ^ The calculated variance (`Nothing` implies it is
                  -- calculated from the supplied map)
                  -> M.Map Double Rational
                  -- ^ The supplied CDF
                  -> EmpiricalCDF
makeEmpiricalCDF' nSamples mMean mVar m
  = EmpiricalCDF
    { _ecdf        = \case
        x | M.null m  -> error "makeEmpiricalCDF: no tangible mass"
          | x <= lwb  -> 0
          | x >= upb  -> mass
          | otherwise -> snd . fromJust $ M.lookupLE x m
    , _ecdfInverse  = \case
        x | M.null i  -> error "makeEmpiricalCDF: no tangible mass"
          | x <  0    -> error "makeEmpiricalCDF: negative probability"
          | x >  1    -> error "makeEmpiricalCDF: > unit probability"
          | x >  mass -> Nothing
          | otherwise -> fmap snd $  M.lookupLE x i
    , _ecdfPdf      = \case
        x | M.null m  -> error "makeEmpiricalCDF: no tangible mass"
          | x <  lwb  -> ((negate infinity, lwb     ), 0)
          | x >= upb  -> ((upb,             infinity), 0)
          | otherwise -> asPDF m x
    , _ecdfMass     = mass
    , _ecdfSamples  = maybe (M.size m) id nSamples
    , _ecdfMin      = lwb
    , _ecdfMax      = upb
    , _ecdfMean     = maybe mean' id mMean
    , _ecdfVariance = maybe var'  id mVar
    }
  where
    mass = M.findWithDefault 0 upb m
    lwb  = maybe infinity fst $ M.lookupMin m
    upb  = maybe (negate infinity) fst $ M.lookupMax m

    -- construct the inverse iCDF defined over the range of the tangible mass
    -- construct the PDF (closure)
    i = M.fromList . map (\(x,y) -> (y,x)) . M.toAscList $ m

    -- Use the probability mass for each point to weight the value. Note that
    -- needs to be scaled by the intangible mass (if any is present).
    mean' :: Double
    mean' = let
      -- need to scale the probabilities so that the weighted sums represent unit
      -- prob mass.
      sf | mass == 0 = 1 -- don't want exception, no non-zero entries anyway
         | otherwise = recip mass
      -- Construct the pointwise PDF
      df _ []         = []
      df n ((a,b):xs) = (a, b - n) : df b xs
      in
        sum [a * fromRational (b * sf) | (a,b) <- df 0 $ M.toAscList  m]

    var'  = nan

-- The result of the lookup is half-open interval [a,b) over which this is the
-- probability density.
    asPDF :: M.Map Double Rational -> Double -> ((Double, Double), Rational)
    asPDF = (f .) . flip M.lookupLE . M.fromAscList . g .  M.toAscList
      where
        f (Just (l, (u, pd))) = ((l, u), pd)
        f Nothing = error "makeEmpiricalCDF: asPDF impossible!"

        g :: [(Double, Rational)] -> [(Double, (Double, Rational))]
        g [] = []
        g (a:as)
          = g' a as
          where
            -- skip over CDF values that are too close to the initial value
            g' _ [] = []
            g' b@(bk,bv) cs@((ck,cv):cs')
              | (ck  < bk * h')  && (not $ null cs')
                = g' b cs'
              | otherwise
                -- d is probablity density over the interval
                = let d = (cv - bv) / toRational (ck - bk)
                  in d `seq` (bk, (ck, d)) : g cs
        -- The square root of the machine epsilon. Used to mitigate potential
        -- numerical instability issues (see
        -- https://en.wikipedia.org/wiki/Numerical_differentiation#Step_size).
        -- Used to combine successive samples where they are too close for
        -- numerical comfort.
        h' :: Double
        h' = 1 + sqrt epsilon
