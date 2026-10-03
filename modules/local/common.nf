// Shared helpers of the KANIT processes. Every step writes into the shared results
// tree (config/config.yaml paths); a task emits receipt.json, which the tasks that
// depend on it stage as input, so that -resume reruns them when it reruns.

def py(String script) {
    "${params.python} ${projectDir}/scripts/${script}"
}

def receipt(String step, key) {
    "${py('receipt.py')} --step ${step} --key '${key}' > receipt.json"
}

def stubReceipt(String step, key) {
    "${py('receipt.py')} --step ${step} --key '${key}' --stub > receipt.json"
}

def meta(String organism, String antibiotic) {
    [id: "${organism}__${antibiotic}", organism: organism, antibiotic: antibiotic]
}

// config/config.yaml, read for the sizes of the fan-outs (the steps read it themselves)
def kanitConfig() {
    new org.yaml.snakeyaml.Yaml().load(new File(params.config_yaml).text)
}
