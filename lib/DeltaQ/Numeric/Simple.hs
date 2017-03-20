{-# Language FlexibleContexts #-}
{-# LANGUAGE FlexibleInstances #-}
module Simple
  ( DQDelay(..)
  , DQEvaluator(..)
  , ProbabilityMass
  , RationalProbabilityMass
  , Weight
  , Delay
  , GenIO
  , createSystemRandom
  , tangibleProbMass
  , tangibleRandomness
  , minDQ
  , runSampleDQ
  , runSampleDQ'
  , singleSampleDQ
  )
where

import Control.Monad.Primitive
import Control.Monad.Reader
import Data.Maybe
import Statistics.Distribution
import Statistics.Distribution.Uniform
import System.IO.Unsafe
import System.Random.MWC
import System.Random.MWC.Distributions

type ProbabilityMass = Double
type RationalProbabilityMass = Rational
type Weight = Rational
type Delay  = Double

data DQDelay
  = DQD_Uniform0 Delay
  | DQD_Dirac    Delay
  | DQD_Convolve DQDelay DQDelay
  | DQD_Choice   (Weight, DQDelay) (Weight, DQDelay)
  | DQD_Perfection
  | DQD_Bottom
  deriving Show

tangibleProbMass :: DQDelay -> RationalProbabilityMass
tangibleProbMass DQD_Bottom = 0
tangibleProbMass (DQD_Convolve a b)
  = tangibleProbMass a * tangibleProbMass b
tangibleProbMass (DQD_Choice (w0,d0) (w1,d1))
  = let {tw = w0 + w1; p0 = w0 / tw; p1 = w1 / tw}
    in (p0 * tangibleProbMass d0) + (p1 * tangibleProbMass d1)
tangibleProbMass _ = 1

tangibleRandomness :: DQDelay -> Maybe DQDelay
tangibleRandomness (DQD_Convolve a b)
  | isJust a' && isJust b' = Just $ DQD_Convolve (fromJust a') (fromJust b')
  | isJust a' = a'
  | isJust b' = b'
  | otherwise = Nothing
  where
    a' = tangibleRandomness a
    b' = tangibleRandomness b
tangibleRandomness (DQD_Choice (w0,a) (w1,b))
  | w0 > 0 && w1 > 0 && isJust a' && isJust b'
   = Just $ DQD_Choice (w0, fromJust a') (w1, fromJust b')
  | w0 > 0 && isJust a' = a'
  | w1 > 0 && isJust b' = b'
  | otherwise = Nothing
  where
    a' = tangibleRandomness a
    b' = tangibleRandomness b
tangibleRandomness DQD_Bottom = Nothing
tangibleRandomness x = Just x

minDQ :: DQDelay -> Delay
minDQ (DQD_Uniform0 _)
  = 0
minDQ (DQD_Dirac x)
  = x
minDQ (DQD_Convolve a b)
  = minDQ a + minDQ b
minDQ (DQD_Choice (w0,a) (w1, b))
  | w0 > 0 && w1 > 0
    = (minDQ a) `min` (minDQ b)
  | w0 > 0
    = minDQ a
  | w1 > 0
    = minDQ b
  | otherwise
    = error "minDQ: zero weights everywhere!"
minDQ DQD_Perfection
  = 0
minDQ DQD_Bottom
  = error "minDQ: only works on tangibleRandomness"

class (MonadIO m) => DQEvaluator m where
  getGenerator :: m GenIO

instance DQEvaluator (ReaderT (Gen RealWorld) IO) where
  getGenerator = ask

runSampleDQ :: (MonadIO m) => DQDelay -> m [Maybe Delay]
runSampleDQ d
  = liftIO $ createSystemRandom >>= runReaderT (runSampleDQ' d)

runSampleDQ' :: (DQEvaluator m) => DQDelay -> m [Maybe Delay]
runSampleDQ' d = do
  gen <- getGenerator
  liftIO $ go gen
  where
    go x = do
      v  <- runReaderT (singleSampleDQ d) x
      vs <- unsafeInterleaveIO . go $ x
      return (v : vs)

getSample :: (ContGen d, DQEvaluator m)
          => d
          -> m Delay
getSample d =
  getGenerator >>= liftIO . genContVar d

getBernoulli :: (DQEvaluator m)
             => RationalProbabilityMass
             -> m Bool
getBernoulli p
  = getGenerator >>= liftIO . bernoulli (fromRational p)


singleSampleDQ :: (DQEvaluator m)
               => DQDelay
               -> m (Maybe Delay)
singleSampleDQ (DQD_Uniform0 x)
  = fmap Just $ getSample (uniformDistr 0 x)
singleSampleDQ (DQD_Dirac x)
  = return $ Just x
singleSampleDQ (DQD_Convolve a b) = do
  a' <- singleSampleDQ a
  b' <- singleSampleDQ b
  case (a', b') of
    (Just a'', Just b'')
      -> return $ Just (a'' + b'')
    _ -> return Nothing
singleSampleDQ (DQD_Choice a b) = do
  weightedChoice a b >>= singleSampleDQ
  where
    weightedChoice (w0,x) (w1,y) = do
      t <- getBernoulli (w0 / (w0 + w1))
      return $ if t then x else y
singleSampleDQ DQD_Perfection
  = return $ Just 0
singleSampleDQ DQD_Bottom
  = return Nothing
