package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object SiestaParserSpec extends Specification {
  "SiestaParserTest" >> {
    "test with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/H2O/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/H2O/out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
