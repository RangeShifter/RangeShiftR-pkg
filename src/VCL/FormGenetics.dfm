object frmGenetics: TfrmGenetics
  Left = 0
  Top = 0
  Caption = 'Genetics Parameters'
  ClientHeight = 460
  ClientWidth = 521
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  PixelsPerInch = 96
  TextHeight = 13
  object LabelChromosome: TLabel
    Left = 100
    Top = 173
    Width = 169
    Height = 21
    Margins.Left = 8
    Margins.Right = 8
    Alignment = taCenter
    AutoSize = False
    Caption = ' Chromosome parameters: '
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold, fsItalic]
    ParentFont = False
  end
  object nTraits: TLabeledEdit
    Left = 154
    Top = 24
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 125
    EditLabel.Height = 16
    EditLabel.Caption = 'No. of variable traits: '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Enabled = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 0
  end
  object CBneutral: TCheckBox
    Left = 272
    Top = 28
    Width = 120
    Height = 17
    Caption = 'Neutral genetics'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentFont = False
    TabOrder = 1
    OnClick = CBneutralClick
  end
  object RGploidy: TRadioGroup
    Left = 130
    Top = 75
    Width = 81
    Height = 69
    Caption = 'Genome'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'Diploid'
      'Haploid')
    ParentFont = False
    TabOrder = 2
  end
  object RGarchitecture: TRadioGroup
    Left = 272
    Top = 75
    Width = 185
    Height = 81
    Caption = 'Architecture'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ItemIndex = 0
    Items.Strings = (
      'One chromosome per trait '
      'Read from file')
    ParentFont = False
    TabOrder = 3
    OnClick = RGarchitectureClick
  end
  object edtSize: TLabeledEdit
    Left = 212
    Top = 200
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 167
    EditLabel.Height = 16
    EditLabel.Caption = 'No. of loci per chromosome: '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 4
    Text = '0'
  end
  object edtMutnProb: TLabeledEdit
    Left = 212
    Top = 240
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 121
    EditLabel.Height = 16
    EditLabel.Caption = 'Mutation probability: '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 5
    Text = '0.0001'
  end
  object edtXoverProb: TLabeledEdit
    Left = 212
    Top = 280
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 129
    EditLabel.Height = 16
    EditLabel.Caption = 'Crossover probability: '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 6
    Text = '0.5'
  end
  object edtAlleleSD: TLabeledEdit
    Left = 212
    Top = 320
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 99
    EditLabel.Height = 16
    EditLabel.Caption = 'Initial allele s.d.: '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 7
    Text = '0.1'
  end
  object edtMutnSD: TLabeledEdit
    Left = 212
    Top = 360
    Width = 57
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 83
    EditLabel.Height = 16
    EditLabel.Caption = 'Mutation s.d.: '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = []
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 8
    Text = '0.1'
  end
  object BtnOK: TButton
    Left = 336
    Top = 410
    Width = 75
    Height = 25
    Caption = 'OK'
    Default = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ModalResult = 1
    ParentFont = False
    TabOrder = 9
    TabStop = False
    OnClick = BtnOKClick
  end
  object BtnReadFile: TButton
    Left = 336
    Top = 200
    Width = 75
    Height = 25
    Caption = 'Read File'
    Default = True
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ModalResult = 1
    ParentFont = False
    TabOrder = 10
    OnClick = BtnReadFileClick
  end
  object OpenDialog1: TOpenDialog
    Left = 416
    Top = 296
  end
end
