package eu.nomad_lab.parsers

import org.specs2.mutable.Specification

object SiestaParserSpec extends Specification {
  "SiestaParserTest" >> {
    "test H2O with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/H2O/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test H2O with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/H2O/out", "json") must_== ParseResult.ParseSuccess
    }
    "test H2O-relax with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/H2O-relax/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test H2O-relax with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/H2O-relax/out", "json") must_== ParseResult.ParseSuccess
    }
    "test Al-bulk with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/Al-bulk/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test Al-bulk with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/Al-bulk/out", "json") must_== ParseResult.ParseSuccess
    }
    "test Al-slab with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/Al-slab/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test Al-slab with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/Al-slab/out", "json") must_== ParseResult.ParseSuccess
    }
    "test Fe with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/Fe/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test Fe with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/Fe/out", "json") must_== ParseResult.ParseSuccess
    }
    "test MgO with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/MgO/out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test MgO with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/MgO/out", "json") must_== ParseResult.ParseSuccess
    }
    "test smeagol-Au-leads with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/smeagol-Au-leads/Au.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test smeagol-Au-leads with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/smeagol-Au-leads/Au.out", "json") must_== ParseResult.ParseSuccess
    }
    "test smeagol-Au-scregion with json-events" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/smeagol-Au-scregion/Au.out", "json-events") must_== ParseResult.ParseSuccess
    }
    "test smeagol-Au-scregion with json" >> {
      ParserRun.parse(SiestaParser, "parsers/siesta/test/examples/smeagol-Au-scregion/Au.out", "json") must_== ParseResult.ParseSuccess
    }
  }
}
