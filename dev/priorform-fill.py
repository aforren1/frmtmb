# Paste the generated blocks of the lane's logs into
# dev/priorform-findings.md verbatim, so no count in it is typed.
#   python3 dev/priorform-fill.py
import re

doc_path = 'dev/priorform-findings.md'
doc = open(doc_path, encoding='utf-8').read()


def block(log, tag):
    text = open(log, encoding='utf-8').read()
    m = re.search(r'<!-- ' + re.escape(tag) + r':begin -->.*?<!-- '
                  + re.escape(tag) + r':end -->', text, re.S)
    if not m:
        raise SystemExit('no block ' + tag + ' in ' + log)
    return m.group(0)


def line(log, prefix):
    for ln in open(log, encoding='utf-8').read().splitlines():
        if ln.startswith(prefix):
            return ln
    raise SystemExit('no line ' + prefix + ' in ' + log)


def replace_block(doc, key, new):
    pat = re.compile(r'<!-- ' + re.escape(key) + r':begin -->.*?<!-- '
                     + re.escape(key) + r':end -->', re.S)
    if pat.search(doc):
        return pat.sub(lambda m: new, doc)
    return doc.replace(key.upper().replace('PRIORFORM-', ''), new)


doc = replace_block(doc, 'priorform-ledger-ref',
                    block('dev/priorform-ledger-ref-log.txt',
                          'priorform-ledger-ref'))
doc = doc.replace('LEDGER_REF', block('dev/priorform-ledger-ref-log.txt',
                                      'priorform-ledger-ref'))
doc = replace_block(doc, 'priorform-ledger-lane',
                    block('dev/priorform-ledger-lane-log.txt',
                          'priorform-ledger-lane'))
doc = doc.replace('LEDGER_LANE', block('dev/priorform-ledger-lane-log.txt',
                                       'priorform-ledger-lane'))
seefail = ('<!-- priorform-seefail:begin -->\n'
           'base build: ' + line('dev/priorform-seefail-ref-log.txt',
                                 'PRIORFORM ') + '\n\n'
           'lane build: ' + line('dev/priorform-seefail-lane-log.txt',
                                 'PRIORFORM ') + '\n'
           '<!-- priorform-seefail:end -->')
doc = replace_block(doc, 'priorform-seefail', seefail)
doc = doc.replace('SEEFAIL', seefail)
doc = replace_block(doc, 'priorform-singular',
                    block('dev/priorform-singular-log.txt',
                          'priorform-singular'))
doc = replace_block(doc, 'priorform-falsealarm',
                    block('dev/priorform-falsealarm-log.txt',
                          'priorform-falsealarm'))

suites = []
logs = [('frmtmb (core)', 'dev/priorform-core-suite-log.txt')] + [
    (e, 'dev/priorform-ext-' + e + '-log.txt') for e in
    ['frmtmb.coupling', 'frmtmb.eam', 'frmtmb.latent', 'frmtmb.learn',
     'frmtmb.ode', 'frmtmb.sample', 'frmtmb.spline']]
try:
    for name, log in logs:
        text = open(log, encoding='utf-8').read()
        m = re.search(r'---- GENERATED COUNTS.*?---- END GENERATED COUNTS ----',
                      text, re.S)
        if not m:
            raise SystemExit('no counts in ' + log)
        suites.append(name + ', ' + log + ':\n\n```\n' + m.group(0) + '\n```')
    s_block = ('<!-- priorform-suites:begin -->\n' + '\n\n'.join(suites)
               + '\n<!-- priorform-suites:end -->')
    doc = replace_block(doc, 'priorform-suites', s_block)
    doc = doc.replace('\nSUITES\n', '\n' + s_block + '\n')
except FileNotFoundError as e:
    print('suites pending:', e)

doc = replace_block(doc, 'priorform-punch-order',
                    block('dev/priorform-punch-order-log.txt',
                          'priorform-punch-order'))
rv = open('dev/priorform-punch-revert-log.txt', encoding='utf-8').read()
rv_block = ('<!-- priorform-punch-revert:begin -->\n```\n' + rv.strip()
            + '\n```\n<!-- priorform-punch-revert:end -->')
doc = replace_block(doc, 'priorform-punch-revert', rv_block)

rv2 = open('dev/priorform-punch2-revert-log.txt', encoding='utf-8').read()
doc = replace_block(doc, 'priorform-punch2-revert',
                    '<!-- priorform-punch2-revert:begin -->\n```\n'
                    + rv2.strip()
                    + '\n```\n<!-- priorform-punch2-revert:end -->')

open(doc_path, 'w', encoding='utf-8', newline='\n').write(doc)
print('filled')
