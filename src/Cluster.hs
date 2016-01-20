{- Cluster
by Gregory Schwartz

-- | Collects the functions pertaining to the clustering of sequences by hamming
distance identity
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ViewPatterns #-}

module Cluster
    ( clusterIdentity
    , groupBy'
    ) where

-- Standard
import qualified Data.Map as Map
import qualified Data.Sequence as Seq
import qualified Data.Foldable as F

-- Cabal
import qualified Data.Text as T
import Data.Fasta.Text

-- Local
import Types

-- | Like hamming, but similarity rather than distance
negHamming :: T.Text -> T.Text -> Double
negHamming x = sum . map (\(x, y) -> if x == y then 1 else 0) . T.zip x

-- | Takes in an identity and fasta sequences and returns the sequences
-- grouped together by hamming distance identity
clusterIdentity :: Identity
                -> Seq.Seq FastaSequence
                -> Map.Map Int (Seq.Seq FastaSequence)
clusterIdentity identity = Map.fromList
                         . zip [1..]
                         . clusterIdentityGo identity

-- | Keep comparing clusters until no more fusions (no change in size) make
-- sense
clusterIdentityGo :: Identity
                  -> Seq.Seq FastaSequence
                  -> [Seq.Seq FastaSequence]
clusterIdentityGo identity = F.toList
                           . groupBy' (compareSeqs identity)

-- | Group together by all pairings rather than adjacent. Altered from
-- lyxia's original to use sequences.
groupBy' :: (a -> a -> Bool) -> Seq.Seq a -> Seq.Seq (Seq.Seq a)
groupBy' _ (Seq.null -> True) = Seq.empty
groupBy' f (Seq.viewl -> x Seq.:< xs) = eqX Seq.<| groupBy' f neqX
  where
    (!eqX, !neqX) = eqTo f Seq.empty (Seq.singleton x) xs

-- | Helper to groupBy'. Altered from lyxia's original to use sequences.
eqTo :: (a -> a -> Bool)
     -> Seq.Seq a
     -> Seq.Seq a
     -> Seq.Seq a
     -> (Seq.Seq a, Seq.Seq a)
eqTo _ acc (Seq.null -> True) zs = (acc, zs)
eqTo f acc (Seq.viewl -> x Seq.:< xs) zs =
    eqTo f (x Seq.<| acc) (eqX Seq.>< xs) neqX
  where
    (!eqX, !neqX) = Seq.partition (f x) zs

-- | Either the sequences are similar or not
compareSeqs :: Identity -> FastaSequence -> FastaSequence -> Bool
compareSeqs identity x y = getIdentity x y >= identity

-- | Get the identity between two sequences
getIdentity :: FastaSequence -> FastaSequence -> Identity
getIdentity xs ys = Identity
                  . (* 100)
                  . (/ (fromIntegral . T.length . fastaSeq $ xs))
                  . negHamming (fastaSeq xs)
                  $ (fastaSeq ys)
