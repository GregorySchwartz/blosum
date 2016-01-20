{- Main
by Gregory Schwartz

-- | Computes a BLOSUM matrix, but takes into account alternate amino acids
-}
module Main where

-- Standard
import Data.Maybe
import qualified System.IO as IO
import qualified Data.Sequence as Seq

-- Cabal
import Options.Applicative
import Data.Fasta.Text
import Pipes
import qualified Pipes.Prelude as P
import qualified Pipes.Text.IO as PT

-- Local
import Types
import Cluster

-- Command line arguments
data Options = Options { input    :: Maybe String
                       , identity :: Double
                       , output   :: Maybe String
                       }

-- Command line options
options :: Parser Options
options = Options
      <$> optional ( strOption
          ( long "input"
         <> short 'i'
         <> metavar "FILE"
         <> help "The input fasta file with blocks labeled in a field"
          )
        )
      <*> option auto
          ( long "identity"
         <> short 'd'
         <> metavar "[100] | PERCENT"
         <> value 100
         <> help "The percent identity of two sequences for clustering.\
                 \ BLOSUM's most famous identity is 62" )
      <*> optional ( strOption
          ( long "output"
         <> short 'o'
         <> metavar "FILE"
         <> help "The output file containing the matrix"
          )
        )

blossum :: Options -> IO ()
blossum opts = do
    h <- case input opts of
            Nothing  -> return IO.stdin
            (Just x) -> IO.openFile x IO.ReadMode
    fastaList <- runEffect $ P.toListM $ pipesFasta . PT.fromHandle $ h

    print $ clusterIdentity (Identity $ identity opts) . Seq.fromList $ fastaList

main :: IO ()
main = execParser opts >>= blossum
  where
    opts = info (helper <*> options)
      ( fullDesc
     <> progDesc "Computes a BLOSUM matrix, but takes into account\
                 \ alternate amino acids"
      )
