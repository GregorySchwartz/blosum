{- Matrix
by Gregory Schwartz

-- | Collects the functions pertaining to the creation of the BLOSSUM
-- matrix
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ViewPatterns #-}

module Matrix
    ( removeGaps
    , getClusterFrequencyMap
    , getBlockMap
    , joinBlockMaps
    , getBlossum
    ) where

-- Standard
import Data.Maybe
import Data.Tuple
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import qualified Data.Foldable as F
import Control.Applicative
import Control.Lens

-- Cabal
import Data.Fasta.Text

-- Local
import Types
import Utility

-- | Zip positions into fasta sequences with a certain size (number of
-- sequences in a cluster) of AA pairs
zipSize :: Frequency -> [AA] -> [(AA, Frequency)]
zipSize x = flip zip [x,x..]

-- | Zip positions into sequences to get position pairs
zipPosition :: [a] -> [(Position, a)]
zipPosition = zip [1..]

-- | Get the frequencies of amino acid pairs for each position in a cluster
getClusterFrequencyMap :: Seq.Seq FastaSequence -> ClusterFrequencyMap
getClusterFrequencyMap xs = ClusterFrequencyMap
                          . Map.map summarize
                          . Map.fromListWith (Seq.><)
                          . concat
                          . F.toList
                          . fmap positionFrequencies
                          $ xs
  where
    summarize           = Seq.fromList . zipSize clusterSize . F.toList
    positionFrequencies = zipPosition
                        . map Seq.singleton
                        . getSeq
    clusterSize         = Frequency . fromIntegral . Seq.length $ xs

-- | Filter gaps out of the map. If no gaps are wanted, remove the entire
-- position
removeGaps :: Bool -> AAMap -> AAMap
removeGaps gapFlag = AAMap
                   . filterGapKey
                   . Map.map filterGapKey
                   . removeAll gapFlag
                   . unAAMap
  where
    removeAll False x = x
    removeAll True x  =
        if or
         . map
           ( \gap -> Map.member gap x
                  || (or . Map.elems . Map.map (Map.member gap) $ x)
           )
         $ gaps
            then Map.empty
            else x
    filterGapKey    = Map.filterWithKey (\k _ -> notElem k gaps)
    gaps            = map AA "-."

collectPairs :: (Ord a, Num b, Fractional b) => Seq.Seq (Seq.Seq (a, b))
                                             -> Seq.Seq (a, (a, b))
                                             -> Seq.Seq (a, (a, b))
collectPairs (Seq.null -> True) !ys         = ys
collectPairs (Seq.viewl -> x Seq.:< xs) !ys =
    collectPairs xs $ comparisons Seq.>< flippedComparisons Seq.>< ys
  where
    flippedComparisons = Seq.filter (\a -> fst a /= (fst . snd $ a))
                       . fmap (\ (!a, (!b, !c)) -> (b, (a, c)))
                       $ comparisons
    comparisons        = F.asum . fmap (pairs x) $ xs
    pairs as bs        = (\ (!a, !b) (!c, !d)
                         -> (a, (c, (1 / (b * d))))
                         )
                     <$> as
                     <*> bs

-- Convert a sequence of clusters of sequences of AAs into an amino acid
-- map. We do such a convoluted method in order to make sure we aren't
-- comparing within a cluster, only between clusters.
toAAMap :: Seq.Seq (Seq.Seq (AA, Frequency)) -> AAMap
toAAMap = AAMap
        . Map.fromListWith (Map.unionWith (+))
        . F.toList
        . fmap (over _2 (uncurry Map.singleton))
        . flip collectPairs Seq.empty

-- | Get the frequency matrix from a list of frequency maps from clusters.
-- We no longer care about positions after this.
getBlockMap :: Bool -> [ClusterFrequencyMap] -> BlockMap
getBlockMap gapFlag = BlockMap
                    . mconcat
                    . Map.elems
                    . Map.map (removeGaps gapFlag . toAAMap)
                    . Map.unionsWith (Seq.><)
                    . map (Map.map Seq.singleton) -- To separate different clusters
                    . map unClusterFrequencyMap

-- | Join together all frequency maps into a single frequency map.
joinBlockMaps :: [BlockMap] -> FrequencyMap
joinBlockMaps = FrequencyMap . mconcat . map unBlockMap

-- | Get the blossum matrix of each AA entry
getBlossum :: FrequencyMap -> Blossum
getBlossum (FrequencyMap (AAMap frequencyMap)) =
    Blossum
        . Map.mapWithKey (\k -> Map.mapWithKey (\l _ -> blossum k l))
        $ frequencyMap
  where
    blossum x y = BlossumVal . round $ 2 * logBase 2 (q x y / e x y)
    e x y       = if x == y then p x * p y else 2 * p x * p y
    p x         = q x x
                + (sum (map (q x) . filter (/= x) . Map.keys $ frequencyMap) / 2)
    q x y       = (\(Frequency a) -> a)
                $ (lookZero y (lookMap x $ frequencyMap))
                / qDenom x y
    lookZero k = fromMaybe 0 . Map.lookup k
    lookMap k  = fromMaybe Map.empty . Map.lookup k
    qDenom x y
        | x == y    = sumMap . lookMap x $ frequencyMap
        | otherwise = (sumMap . lookMap x $ frequencyMap)
                    + ( sumMap
                      . Map.filterWithKey (\k _ -> k /= x)
                      . lookMap y
                      $ frequencyMap
                      )
