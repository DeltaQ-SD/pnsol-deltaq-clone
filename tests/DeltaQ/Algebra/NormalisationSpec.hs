{-# LANGUAGE DerivingStrategies #-}
{-# LANGUAGE FlexibleInstances #-}
{-# LANGUAGE GeneralizedNewtypeDeriving #-}

module DeltaQ.Algebra.NormalisationSpec where

import Control.Monad (liftM)
import Debug.Trace (traceShow)
import DeltaQ.Algebra (DeltaQ (..), normaliseDeltaQ)
import DeltaQ.Algebra.DelayModel.SimpleUniform (SimpleUniform (..))
import Test.Hspec (Spec)
import Test.Hspec.QuickCheck (prop)
import Test.QuickCheck (Arbitrary (..), Gen, Property, choose, forAll, oneof, suchThat)

spec :: Spec
spec = do
    prop "a [0⇋1] b ⇒ b" prop_p_rule_zero_choice
    prop "a [1⇋0] b ⇒ a" prop_p_rule_one_choice
    prop "⊥ [p⇋q] ⊥ ⇒ ⊥, ∀p,q" prop_p_rule_bottom
    prop "∅ [p⇋q] ∅ ⇒ ∅, ∀p,q" prop_p_rule_unit
    prop "x [p⇋q] ⊥ ⇒ ⊥ [q⇋p] x" prop_p_norm_bottom
    prop "∅ [p⇋q] x ⇒ x [q⇋p] ∅" prop_p_norm_unit
    prop "⊥ ⊕ x ⇒ ⊥, ∀x" prop_c_bottom_a
    prop "x ⊕ ⊥ ⇒ ⊥, ∀x" prop_c_bottom_b
    prop "∅ ⊕ x ⇒ x, ∀x" prop_c_unit_a
    prop "x ⊕ ∅ ⇒ x, ∀x" prop_c_unit_b

-- testGroup "⊥ promotion"
-- test "(⊥ [p⇋q] x) ⊕ y ⇒ ⊥ [p⇋q] (x ⊕ y)" prop_bottom_promotion_a
-- , test "x ⊕ (⊥ [p⇋q] y) ⇒ ⊥ [p⇋q] (x ⊕ y)" prop_bottom_promotion_b

-- testGroup "∅ elimination"
-- test "(x [p⇋q] ∅) ⊕ y ⇒ (x ⊕ y) [p⇋q] y" prop_unit_elimination_a
-- , test "x ⊕ (y [p⇋q] ∅) ⇒ (x ⊕ y) [p⇋q] x" prop_unit_elimination_b

-- testGroup "⊥ concatentation"
-- test "⊥ [p⇋q] (⊥ [r⇋s] x) ⇒ ⊥ [r*(p+q)+s*p⇋(p+q)*(r+s)] x"
--               prop_bottom_concatenation

-- testGroup "∅ demotion"
-- test "(x [p⇋q] ∅) [r⇋s] y ⇒ x [p*r⇋]

prop_p_rule_zero_choice :: RationalDeltaQ -> Property
prop_p_rule_zero_choice _ =
    forAll tc $ \x ->
        case x of
            (ProbChoice _ _ b) -> (normaliseDeltaQ x == normaliseDeltaQ b)
            _ -> error $ "prop_p_rule_zero_choice generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = ProbChoice 0 <$> arbitrary <*> arbitrary

prop_p_rule_one_choice :: RationalDeltaQ -> Property
prop_p_rule_one_choice _ =
    forAll tc $ \x ->
        case x of
            (ProbChoice _ a _) -> (normaliseDeltaQ x == normaliseDeltaQ a)
            _ -> error $ "prop_p_rule_one_choice generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = ProbChoice 1 <$> arbitrary <*> arbitrary

prop_p_rule_bottom :: RationalDeltaQ -> Property
prop_p_rule_bottom _ =
    forAll tc $ \x ->
        case x of
            (ProbChoice _ Bottom Bottom) -> (normaliseDeltaQ x == Bottom)
            _ -> error $ "prop_p_rule_bottom generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = liftM (\y -> ProbChoice y Bottom Bottom) arbitrary

prop_p_rule_unit :: RationalDeltaQ -> Property
prop_p_rule_unit _ =
    forAll tc $ \x ->
        case x of
            (ProbChoice _ Unit Unit) -> (normaliseDeltaQ x == Unit)
            _ -> error $ "prop_p_rule_unit generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = liftM (\y -> ProbChoice y Unit Unit) arbitrary

prop_p_norm_bottom :: RationalDeltaQ -> Property
prop_p_norm_bottom _ =
    forAll tc $ \z ->
        case z of
            (ProbChoice _p _x Bottom) -> case normaliseDeltaQ z of
                (ProbChoice _ Bottom _y) -> True
                Bottom -> normaliseDeltaQ _x == Bottom
                _ -> False
            _ -> error $ "prop_p_norm_bottom generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = (\p x -> ProbChoice p x Bottom) <$> not0or1 <*> arbitrary
    not0or1 = arbitrary `suchThat` (\p -> p > 0 && p < 1)

prop_p_norm_unit :: RationalDeltaQ -> Property
prop_p_norm_unit _ =
    forAll tc $ \z ->
        case z of
            (ProbChoice _p Unit _y) ->
                case normaliseDeltaQ z of
                    Unit -> normaliseDeltaQ _y == Unit
                    (ProbChoice _ _ y) -> endsInUnit $ normaliseDeltaQ y
                    _z -> traceShow ("prob_p_norm_unit", z, _z) False
            _ -> error $ "prop_p_norm_unit generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = (\p y -> ProbChoice p Unit y) <$> not0or1 <*> arbitrary
    not0or1 = arbitrary `suchThat` (\p -> p > 0 && p < 1)

    endsInUnit Unit = True
    endsInUnit (ProbChoice _ _ y) = endsInUnit y
    endsInUnit _ = False

prop_c_bottom_a :: RationalDeltaQ -> Property
prop_c_bottom_a _ =
    forAll tc $ \z ->
        case z of
            (Convolve Bottom x) ->
                case normaliseDeltaQ z of
                    Bottom -> True
                    _z -> traceShow (name, z, _z) False
            _ -> error $ name ++ " generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = liftM (\x -> Convolve Bottom x) arbitrary
    name = "prop_c_bottom_a"

prop_c_bottom_b :: RationalDeltaQ -> Property
prop_c_bottom_b _ =
    forAll tc $ \z ->
        case z of
            (Convolve x Bottom) ->
                case normaliseDeltaQ z of
                    Bottom -> True
                    _z -> traceShow (name, z, _z) False
            _ -> error $ name ++ " generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = liftM (\x -> Convolve x Bottom) arbitrary
    name = "prop_c_bottom_b"

prop_c_unit_a :: RationalDeltaQ -> Property
prop_c_unit_a _ =
    forAll tc $ \z ->
        case z of
            (Convolve Unit x) ->
                case normaliseDeltaQ z == normaliseDeltaQ x of
                    True -> True
                    False -> traceShow (name, z, normaliseDeltaQ z) False
            _ -> error $ name ++ " generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = liftM (\x -> Convolve Unit x) arbitrary
    name = "prop_c_unit_a"

prop_c_unit_b :: RationalDeltaQ -> Property
prop_c_unit_b _ =
    forAll tc $ \z ->
        case z of
            (Convolve x Unit) ->
                case normaliseDeltaQ z == normaliseDeltaQ x of
                    True -> True
                    False -> traceShow (name, z, normaliseDeltaQ z) False
            _ -> error $ name ++ " generator error"
  where
    tc :: Gen RationalDeltaQ
    tc = liftM (\x -> Convolve x Unit) arbitrary
    name = "prop_c_unit_b"

type RationalDeltaQ = DeltaQ Probability SimpleUniform NonNeg

-- another interative helper
ndq :: RationalDeltaQ -> RationalDeltaQ
ndq = normaliseDeltaQ

instance Arbitrary RationalDeltaQ where
    arbitrary =
        oneof
            [ return Bottom
            , return Unit
            , liftM (Delay . DiracDelta) arbitrary
            , liftM (Delay . UniformD) arbitrary
            , Convolve <$> arbitrary <*> arbitrary
            , ProbChoice <$> arbitrary <*> arbitrary <*> arbitrary
            ]

newtype Probability = Probability Rational
    deriving newtype (Eq, Ord, Show, Num, Fractional, Real)

instance Bounded Probability where
    minBound = Probability 0
    maxBound = Probability 1

newtype NonNeg = NonNeg Rational
    deriving (Eq, Ord, Show, Num, Fractional, Real)
instance Bounded NonNeg where
    minBound = NonNeg 0
    maxBound = NonNeg 1000 -- just an arbitrary choice

instance Arbitrary Probability where
    arbitrary = choose (0, 1) >>= return . Probability . fromRational . (toRational :: Double -> Rational)

instance Arbitrary NonNeg where
    arbitrary = choose (0, 1000) >>= return . NonNeg . fromRational . (toRational :: Double -> Rational)
