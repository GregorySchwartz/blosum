{- Main
by Gregory Schwartz

-- | Computes a BLOSUM matrix, but takes into account alternate amino acids
-}
module Main where

-- Standard
import Data.Maybe
import qualified Data.Map.Strict as Map
import qualified Data.Sequence as Seq
import qualified System.IO as IO
import qualified Data.List.Split as Split
import Debug.Trace

-- Cabal
import qualified Data.Text as T
import qualified Data.Text.IO as T
import Data.Fasta.Text
import Pipes
import qualified Pipes.Prelude as P
import qualified Pipes.Text.IO as PT
import Options.Applicative

-- Local
import Types
import Utility
import Cluster
import Matrix
import Print

-- Command line arguments
data Options = Options { input      :: Maybe String
                       , identity   :: Double
                       , blockField :: Maybe Int
                       , gapFlag    :: Bool
                       , output     :: Maybe String
                       }

-- Command line options
options :: Parser Options
options = Options
      <$> optional ( strOption
          ( long "input"
         <> short 'i'
         <> metavar "FILE"
         <> help "The input fasta file with blocks labeled in a field.\
                 \ If empty, then assumed to use stdin to send the many, many\
                 \ filenames, where each file is a block."
          )
        )
      <*> option auto
          ( long "identity"
         <> short 'd'
         <> metavar "[100] | PERCENT"
         <> value 100
         <> help "The percent identity of two sequences for clustering.\
                 \ BLOSUM's most famous identity is 62" )
      <*> optional ( option auto
          ( long "field"
         <> short 'f'
         <> metavar "FIELD"
         <> help "The field in the fasta headers (1 indexed) which says\
                 \ the block it belongs to. Assumed to already be next to\
                 \ each other, so the block members should not be scattered" )
        )
      <*> switch
          ( long "remove-gaps"
         <> short 'g'
         <> help "Whether to remove positions that have a gap" )
      <*> optional ( strOption
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> help "The output file containing the matrix"
          )
        )

-- | Get a frequency map from a single file
getFrequencyMapSingleFile :: Options -> FilePath -> IO FrequencyMap
getFrequencyMapSingleFile opts file = do
    h         <- IO.openFile file IO.ReadMode
    fastaList <- runEffect $ P.toListM $ pipesFasta . PT.fromHandle $ h

    let clusters      = getClusterIdentity (Identity $ identity opts)
        clusterMaps   = Map.elems
                      . Map.map getClusterFrequencyMap
                      . unClusterMap
        blockMaps     = fmap (getBlockMap (gapFlag opts))
                      . fmap clusterMaps
                      . fmap clusters
                      . fmap Seq.fromList
                      . groupBlocks (fmap Field $ blockField opts)
                      $ fastaList
        frequencyMap  = joinBlockMaps blockMaps

    return frequencyMap

-- | Get a block from a filepath
getBlock :: Options -> FilePath -> IO BlockMap
getBlock opts file = do
    h         <- IO.openFile file IO.ReadMode
    fastaList <- runEffect $ P.toListM $ pipesFasta . PT.fromHandle $ h

    let clusters   = getClusterIdentity (Identity $ identity opts)
        clusterMap = Map.elems
                   . Map.map getClusterFrequencyMap
                   . unClusterMap
        blockMap   = getBlockMap (gapFlag opts)
                   . clusterMap
                   . clusters
                   . Seq.fromList
                   $ fastaList

    return blockMap

blossum :: Options -> IO ()
blossum opts = do

    frequencyMap <-
        case input opts of
            Nothing  -> fmap FrequencyMap
                      . P.fold (<>) (AAMap Map.empty) id
                      $ P.stdinLn
                    >-> P.mapM (fmap unBlockMap . getBlock opts)
            (Just x) -> getFrequencyMapSingleFile opts x

    let blossumMatrix = getBlossum frequencyMap

    -- Save results
    case output opts of
        Nothing  -> T.putStrLn . printBlossum $ blossumMatrix
        (Just x) -> T.writeFile x . printBlossum $ blossumMatrix

main :: IO ()
main = execParser opts >>= blossum
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Computes a BLOSUM matrix, but takes into account\
                 \ alternate amino acids"
      )
