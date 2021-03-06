{- Matrix
by Gregory Schwartz

-- | Collects the functions pertaining to the creation of the BLOSUM
-- matrix
-}

{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE ViewPatterns #-}

module Matrix
    ( removeChar
    , getClusterFrequencyMap
    , getBlockMap
    , joinBlockMaps
    , getBlosum
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
import Data.Fasta.Text (FastaSequence)

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
removeChar :: Bool -> Maybe [AA] -> AAMap -> AAMap
removeChar _ Nothing = id
removeChar allFlag (Just badChars) = AAMap
                                   . filterBadKey
                                   . Map.map filterBadKey
                                   . removeAll allFlag
                                   . unAAMap
  where
    removeAll False x = x
    removeAll True x  =
        if or
         . map
           ( \bad -> Map.member bad x
                  || (or . Map.elems . Map.map (Map.member bad) $ x)
           )
         $ badChars
            then Map.empty
            else x
    filterBadKey    = Map.filterWithKey (\k _ -> notElem k badChars)

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
getBlockMap :: Bool -> Maybe [AA] -> [ClusterFrequencyMap] -> BlockMap
getBlockMap allFlag badChars = BlockMap
                             . mconcat
                             . Map.elems
                             . Map.map (removeChar allFlag badChars . toAAMap)
                             . Map.unionsWith (Seq.><)
                             . map (Map.map Seq.singleton) -- To separate different clusters
                             . map unClusterFrequencyMap

-- | Join together all frequency maps into a single frequency map.
joinBlockMaps :: [BlockMap] -> FrequencyMap
joinBlockMaps = FrequencyMap . mconcat . map unBlockMap

-- | Get the blosum matrix of each AA entry
getBlosum :: FrequencyMap -> Blosum
getBlosum (FrequencyMap (AAMap frequencyMap)) =
    Blosum
        . Map.mapWithKey (\k -> Map.mapWithKey (\l _ -> blosum k l))
        $ frequencyMap
  where
    blosum x y = BlosumVal . round $ 2 * logBase 2 (q x y / e x y)
    e x y      = if x == y then p x * p y else 2 * p x * p y
    p x        = q x x
               + (sum (map (q x) . filter (/= x) . Map.keys $ frequencyMap) / 2)
    q x y      = (\(Frequency a) -> a)
               $ (lookZero y (lookMap x $ frequencyMap))
               / qDenom
    qDenom = (/ 2)
           $ (sumMap . Map.map sumMap $ frequencyMap)
           + ( sumMap
             . Map.mapWithKey
               (\k1 -> sumMap . Map.filterWithKey (\k2 _ -> k1 == k2))
             $ frequencyMap
             )
    lookZero k = fromMaybe 0 . Map.lookup k
    lookMap k  = fromMaybe Map.empty . Map.lookup k
