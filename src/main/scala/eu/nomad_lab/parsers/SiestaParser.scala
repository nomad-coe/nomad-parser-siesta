package eu.nomad_lab.parsers

import eu.{ nomad_lab => lab }
import eu.nomad_lab.DefaultPythonInterpreter
import org.{ json4s => jn }
import scala.collection.breakOut

object SiestaParser extends SimpleExternalParserGenerator(
  name = "SiestaParser",
  parserInfo = jn.JObject(
    ("name" -> jn.JString("SiestaParser")) ::
      ("parserId" -> jn.JString("SiestaParser" + lab.SiestaVersionInfo.version)) ::
      ("versionInfo" -> jn.JObject(
        ("nomadCoreVersion" -> jn.JObject(lab.NomadCoreVersionInfo.toMap.map {
          case (k, v) => k -> jn.JString(v.toString)
        }(breakOut): List[(String, jn.JString)])) ::
          (lab.SiestaVersionInfo.toMap.map {
            case (key, value) =>
              (key -> jn.JString(value.toString))
          }(breakOut): List[(String, jn.JString)])
      )) :: Nil
  ),
  mainFileTypes = Seq("text/.*"),
  //TODO: Update the replacement string (mainFileRe)
  mainFileRe = """(Siesta Version: siesta-|SIESTA [0-9]\.[0-9]\.[0-9])""".r,
  cmd = Seq(DefaultPythonInterpreter.pythonExe(), "${envDir}/parsers/siesta/parser/parser-siesta/main.py",
    "${mainFilePath}"),
  resList = Seq(
    "parser-siesta/main.py",
    "parser-siesta/inputvars.py",
    "parser-siesta/util.py",
    "parser-siesta/setup_paths.py",
    "nomad_meta_info/public.nomadmetainfo.json",
    "nomad_meta_info/common.nomadmetainfo.json",
    "nomad_meta_info/meta_types.nomadmetainfo.json",
    "nomad_meta_info/siesta.autogenerated.nomadmetainfo.json",
    "nomad_meta_info/siesta.nomadmetainfo.json"
  ) ++ DefaultPythonInterpreter.commonFiles(),
  dirMap = Map(
    "parser-siesta" -> "parsers/siesta/parser/parser-siesta",
    "nomad_meta_info" -> "nomad-meta-info/meta_info/nomad_meta_info"
  ) ++ DefaultPythonInterpreter.commonDirMapping()
)
