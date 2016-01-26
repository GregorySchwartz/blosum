{- Matrix
by Gregory Schwartz

-- | Collects the functions pertaining to the creation of the BLOSSUM
-- matrix
-}

{-# LANGUAGE BangPatterns #-}

module Matrix
    ( removeGaps
    , getClusterFrequencyMap
    , getBlockMap
    , joinBlockMaps
    , getBlossum
    ) where

-- Standard
import Data.Maybe
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
                          . Map.map (Seq.fromList . nub' . F.toList)
                          . Map.fromListWith (Seq.><)
                          . concat
                          . F.toList
                          . fmap positionFrequencies
                          $ xs
  where
    positionFrequencies = zipPosition
                        . map Seq.singleton
                        . getSeq
    clusterSize         = Frequency . fromIntegral . Seq.length $ xs

-- | Filter gaps out of the map. If no gaps are wanted, remove the entire
-- position
removeGaps :: Bool -> BlockMap -> BlockMap
removeGaps gapFlag = BlockMap
                   . AAMap
                   . removeAll gapFlag
                   . Map.map filterGapKey
                   . unAAMap
                   . unBlockMap
  where
    removeAll False = id
    removeAll True  = filterGapKey
    filterGapKey    = Map.filterWithKey (\k _ -> notElem k gaps)
    gaps            = map AA "-."

-- | Count all pairs of tuples
collectPairs :: (Ord a) => Seq.Seq a -> Seq.Seq a -> Seq.Seq (a, a)
collectPairs = liftA2 (,)

-- Convert a sequence of clusters of sequences of AAs into an amino acid
-- map. We do such a convoluted method in order to make sure we aren't
-- comparing within a cluster, only between clusters.
toAAMap :: Seq.Seq (Seq.Seq AA) -> AAMap
toAAMap x = AAMap
          . Map.fromListWith (Map.unionWith (+))
          . F.toList
          . fmap (over _2 (flip Map.singleton (Frequency 1)))
          . F.asum
          . cutEnds
          . liftA2 collectPairs x
          $ x

-- | Get the frequency matrix from a list of frequency maps from clusters.
-- We no longer care about positions after this.
getBlockMap :: [ClusterFrequencyMap] -> BlockMap
getBlockMap = BlockMap
            . mconcat
            . Map.elems
            . Map.map toAAMap
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
    q x y       = (\(Frequency x) -> x)
                $ (lookZero y (lookMap x $ frequencyMap))
                / (Map.foldl' (+) 0 $ lookMap x frequencyMap)
    lookZero k = fromMaybe 0 . Map.lookup k
    lookMap k  = fromMaybe Map.empty . Map.lookup k
