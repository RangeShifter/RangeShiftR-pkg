object frmSim: TfrmSim
  Left = 51
  Top = 56
  Width = 777
  Height = 688
  VertScrollBar.ButtonSize = 10
  VertScrollBar.Margin = 5
  VertScrollBar.Tracking = True
  AutoScroll = True
  Caption = 'Simulation Parameters'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'Tahoma'
  Font.Style = []
  OldCreateOrder = False
  Position = poDesigned
  PixelsPerInch = 96
  TextHeight = 13
  object Bevel2: TBevel
    Left = 435
    Top = 10
    Width = 10
    Height = 84
    Shape = bsLeftLine
  end
  object SimulationDesc: TLabel
    Left = 17
    Top = 40
    Width = 186
    Height = 54
    AutoSize = False
    Caption = 
      'This number will define the names of the Output  and Initialisat' +
      'ion files '
    Color = clBtnFace
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -12
    Font.Name = 'Tahoma'
    Font.Style = [fsItalic]
    ParentColor = False
    ParentFont = False
    WordWrap = True
  end
  object CBabsorbing: TCheckBox
    Left = 280
    Top = 70
    Width = 149
    Height = 17
    Alignment = taLeftJustify
    Caption = 'Absorbing boundaries'
    Color = clBtnFace
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    ParentColor = False
    ParentFont = False
    TabOrder = 8
  end
  object PanelOutputs: TPanel
    Left = 10
    Top = 268
    Width = 740
    Height = 344
    BevelInner = bvLowered
    TabOrder = 5
    object RGMap: TRadioGroup
      Left = 10
      Top = 270
      Width = 85
      Height = 72
      Caption = 'Save Maps'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 0
      Items.Strings = (
        'No'
        'Yes')
      ParentFont = False
      TabOrder = 2
      OnClick = RGMapClick
    end
    object PanelMaps: TPanel
      Left = 100
      Top = 270
      Width = 214
      Height = 74
      BevelEdges = []
      BevelOuter = bvNone
      Enabled = False
      TabOrder = 4
      Visible = False
      object edtMap_interval: TLabeledEdit
        Left = 105
        Top = 10
        Width = 40
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 89
        EditLabel.Height = 16
        EditLabel.Caption = 'Every (years):  '
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
        TabOrder = 0
        Text = '10'
      end
      object CBDrawInit: TCheckBox
        Left = 14
        Top = 40
        Width = 211
        Height = 25
        Caption = 'Draw loaded species distribution'
        Color = clBtnFace
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentColor = False
        ParentFont = False
        TabOrder = 1
        WordWrap = True
      end
    end
    object GroupBoxOutputs: TGroupBox
      Left = 10
      Top = 9
      Width = 500
      Height = 255
      Caption = 'Outputs'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 0
      object LabelFreq: TLabel
        Left = 390
        Top = 10
        Width = 59
        Height = 16
        Caption = 'Frequency'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
      end
      object LabelStart: TLabel
        Left = 275
        Top = 10
        Width = 57
        Height = 16
        Caption = 'Start year'
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
      end
      object CBoutRange: TCheckBox
        Left = 16
        Top = 35
        Width = 145
        Height = 17
        Caption = 'Range'
        TabOrder = 0
      end
      object CBoutOcc: TCheckBox
        Left = 16
        Top = 60
        Width = 145
        Height = 17
        Caption = 'Occupancy'
        TabOrder = 21
        OnClick = CBoutOccClick
      end
      object CBoutPop: TCheckBox
        Left = 16
        Top = 85
        Width = 97
        Height = 17
        Caption = 'Populations'
        TabOrder = 1
        OnClick = CBoutPopClick
      end
      object CBoutInd: TCheckBox
        Left = 16
        Top = 110
        Width = 121
        Height = 17
        Caption = 'Individuals'
        TabOrder = 2
        OnClick = CBoutIndClick
      end
      object CBoutGen: TCheckBox
        Left = 16
        Top = 135
        Width = 121
        Height = 17
        Caption = 'Genetics'
        TabOrder = 22
        OnClick = CBoutGenClick
      end
      object CBoutTraitCell: TCheckBox
        Left = 16
        Top = 160
        Width = 145
        Height = 17
        Caption = 'Mean Traits by cells'
        Enabled = False
        TabOrder = 3
        OnClick = CBoutTraitCellClick
      end
      object CBoutTraitRow: TCheckBox
        Left = 16
        Top = 185
        Width = 153
        Height = 17
        Caption = 'Mean Traits by rows'
        Enabled = False
        TabOrder = 4
        OnClick = CBoutTraitRowClick
      end
      object CBoutConnect: TCheckBox
        Left = 16
        Top = 210
        Width = 153
        Height = 17
        Caption = 'Connectivity Matrix'
        Enabled = False
        TabOrder = 5
        OnClick = CBoutConnectClick
      end
      object edtFreqRange: TEdit
        Left = 400
        Top = 30
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 23
        Text = '1'
      end
      object edtFreqOcc: TEdit
        Left = 400
        Top = 55
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 8
        Text = '10'
      end
      object edtFreqPop: TEdit
        Left = 400
        Top = 80
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 9
        Text = '10'
      end
      object edtFreqInd: TEdit
        Left = 400
        Top = 105
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 10
        Text = '10'
      end
      object edtFreqGen: TEdit
        Left = 400
        Top = 130
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 7
        Text = '10'
      end
      object edtFreqTraitCell: TEdit
        Left = 400
        Top = 155
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 11
        Text = '10'
      end
      object edtFreqTraitRow: TEdit
        Left = 400
        Top = 180
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 12
        Text = '10'
      end
      object edtFreqConn: TEdit
        Left = 400
        Top = 205
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 13
        Text = '10'
        Visible = False
      end
      object edtStartRange: TEdit
        Left = 280
        Top = 30
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 14
        Text = '0'
        Visible = False
      end
      object edtStartOcc: TEdit
        Left = 280
        Top = 55
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 15
        Text = '0'
        Visible = False
      end
      object edtStartPop: TEdit
        Left = 280
        Top = 80
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 16
        Text = '0'
      end
      object edtStartInd: TEdit
        Left = 280
        Top = 105
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 17
        Text = '0'
      end
      object edtStartGen: TEdit
        Left = 280
        Top = 130
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 6
        Text = '0'
      end
      object edtStartTraitCell: TEdit
        Left = 280
        Top = 155
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 18
        Text = '0'
      end
      object edtStartTraitRow: TEdit
        Left = 280
        Top = 180
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 19
        Text = '0'
      end
      object edtStartConn: TEdit
        Left = 280
        Top = 205
        Width = 40
        Height = 24
        Alignment = taRightJustify
        Enabled = False
        Font.Charset = DEFAULT_CHARSET
        Font.Color = clWindowText
        Font.Height = -13
        Font.Name = 'Tahoma'
        Font.Style = []
        ParentFont = False
        TabOrder = 20
        Text = '0'
        Visible = False
      end
    end
    object GroupBoxVisual: TGroupBox
      Left = 530
      Top = 9
      Width = 199
      Height = 193
      Caption = 'Dynamic visualisations (slower)'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 1
      object CBVisualLand: TCheckBox
        Left = 8
        Top = 24
        Width = 97
        Height = 17
        Caption = 'Landscape'
        Checked = True
        State = cbChecked
        TabOrder = 0
        OnClick = CBVisualLandClick
      end
      object CBVisualTraits: TCheckBox
        Left = 8
        Top = 116
        Width = 97
        Height = 17
        Caption = 'Mean Traits'
        Enabled = False
        TabOrder = 4
        OnClick = CBVisualTraitsClick
      end
      object CBVisualGraph: TCheckBox
        Left = 8
        Top = 70
        Width = 169
        Height = 17
        Caption = 'Population size (graph) '
        Checked = True
        State = cbChecked
        TabOrder = 2
      end
      object CBVisualPop: TCheckBox
        Left = 8
        Top = 47
        Width = 171
        Height = 17
        Caption = 'Population size (map)'
        Checked = True
        State = cbChecked
        TabOrder = 1
      end
      object CBVisualGrad: TCheckBox
        Left = 8
        Top = 93
        Width = 97
        Height = 17
        Caption = 'Env. Gradient'
        Enabled = False
        TabOrder = 3
      end
      object CBVisualMove: TCheckBox
        Left = 8
        Top = 139
        Width = 161
        Height = 17
        Caption = 'Movement paths'
        Enabled = False
        TabOrder = 5
        OnClick = CBVisualMoveClick
      end
      object edtMovtSlow: TLabeledEdit
        Left = 80
        Top = 160
        Width = 56
        Height = 24
        Alignment = taRightJustify
        EditLabel.Width = 70
        EditLabel.Height = 16
        EditLabel.Caption = 'Slow factor:'
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
        Text = '1'
      end
    end
    object RGTraitsMap: TRadioGroup
      Left = 320
      Top = 270
      Width = 117
      Height = 72
      Caption = 'Save Traits Maps'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 0
      Items.Strings = (
        'No'
        'Yes')
      ParentFont = False
      TabOrder = 3
      Visible = False
      OnClick = RGTraitsMapClick
    end
    object RGGen: TRadioGroup
      Left = 530
      Top = 199
      Width = 199
      Height = 95
      Caption = 'Output genetics for:'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 0
      Items.Strings = (
        'juveniles only'
        'all individuals'
        'adults only')
      ParentFont = False
      TabOrder = 6
      Visible = False
    end
    object edtTraitsMap_Int: TLabeledEdit
      Left = 463
      Top = 315
      Width = 47
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 46
      EditLabel.Height = 32
      EditLabel.BiDiMode = bdLeftToRight
      EditLabel.Caption = 'Every (years):'
      EditLabel.Font.Charset = DEFAULT_CHARSET
      EditLabel.Font.Color = clWindowText
      EditLabel.Font.Height = -13
      EditLabel.Font.Name = 'Tahoma'
      EditLabel.Font.Style = []
      EditLabel.ParentBiDiMode = False
      EditLabel.ParentFont = False
      EditLabel.Transparent = False
      EditLabel.Layout = tlCenter
      EditLabel.WordWrap = True
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 5
      Text = '10'
      Visible = False
    end
    object CBGenCrosstab: TCheckBox
      Left = 530
      Top = 300
      Width = 199
      Height = 17
      Caption = 'Output genetics as cross table'
      Color = clBtnFace
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentColor = False
      ParentFont = False
      TabOrder = 7
      Visible = False
    end
  end
  object BtnOK: TButton
    Left = 489
    Top = 617
    Width = 89
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
    TabOrder = 6
    TabStop = False
    OnClick = BtnOKClick
  end
  object PanelEnvStoch: TPanel
    Left = 10
    Top = 100
    Width = 740
    Height = 162
    BevelInner = bvLowered
    TabOrder = 4
    object Bevel3: TBevel
      Left = 0
      Top = 117
      Width = 740
      Height = 7
      Shape = bsBottomLine
    end
    object LabelAC: TLabel
      Left = 421
      Top = 45
      Width = 92
      Height = 16
      Caption = '0.0 <= ac < 1.0'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      Transparent = False
    end
    object LabelStd: TLabel
      Left = 417
      Top = 95
      Width = 96
      Height = 16
      Caption = '0.0 < std <= 1.0'
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
    end
    object edtAC: TLabeledEdit
      Left = 463
      Top = 20
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 150
      EditLabel.Height = 16
      EditLabel.Caption = 'Temporal autocorrelation:'
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
      TabOrder = 3
      Text = '0.0'
    end
    object edtStd: TLabeledEdit
      Left = 463
      Top = 65
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 195
      EditLabel.Height = 32
      EditLabel.Caption = 'Amplitude (standard deviation for the random normal variable):'
      EditLabel.Font.Charset = DEFAULT_CHARSET
      EditLabel.Font.Color = clWindowText
      EditLabel.Font.Height = -13
      EditLabel.Font.Name = 'Tahoma'
      EditLabel.Font.Style = []
      EditLabel.ParentFont = False
      EditLabel.WordWrap = True
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      LabelPosition = lpLeft
      ParentFont = False
      TabOrder = 4
      Text = '0.25'
    end
    object CBEnvStoch: TCheckBox
      Left = 8
      Top = 12
      Width = 217
      Height = 17
      Caption = 'Environmental Stochasticity'
      Color = clWindowText
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentColor = False
      ParentFont = False
      TabOrder = 0
      OnClick = CBEnvStochClick
    end
    object edtMaxR: TLabeledEdit
      Left = 657
      Top = 50
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 103
      EditLabel.Height = 16
      EditLabel.Caption = 'Max. growth rate:'
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
      TabOrder = 6
      Text = '0.0'
    end
    object edtMinR: TLabeledEdit
      Left = 657
      Top = 20
      Width = 50
      Height = 24
      Alignment = taRightJustify
      EditLabel.Width = 100
      EditLabel.Height = 16
      EditLabel.Caption = 'Min. growth rate:'
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
      TabOrder = 5
      Text = '0.0'
    end
    object CBLocExt: TCheckBox
      Left = 8
      Top = 130
      Width = 217
      Height = 17
      Caption = 'Local extinction probability'
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = [fsBold]
      ParentFont = False
      TabOrder = 7
      OnClick = CBLocExtClick
    end
    object edtLocExt: TEdit
      Left = 217
      Top = 130
      Width = 50
      Height = 24
      Alignment = taRightJustify
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ParentFont = False
      TabOrder = 8
      Text = '0.0'
      Visible = False
    end
    object RGEnvStoch: TRadioGroup
      Left = 8
      Top = 35
      Width = 81
      Height = 71
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 1
      Items.Strings = (
        'Global'
        'Local ')
      ParentFont = False
      TabOrder = 1
      OnClick = RGEnvStochClick
    end
    object RGEnvStochType: TRadioGroup
      Left = 95
      Top = 35
      Width = 146
      Height = 71
      Enabled = False
      Font.Charset = DEFAULT_CHARSET
      Font.Color = clWindowText
      Font.Height = -13
      Font.Name = 'Tahoma'
      Font.Style = []
      ItemIndex = 0
      Items.Strings = (
        'in growth rate'
        'in carrying capacity')
      ParentFont = False
      TabOrder = 2
      OnClick = RGEnvStochTypeClick
    end
  end
  object BtnCancel: TButton
    Left = 584
    Top = 617
    Width = 89
    Height = 25
    Cancel = True
    Caption = 'Cancel'
    ModalResult = 2
    TabOrder = 7
    OnClick = BtnCancelClick
  end
  object edtYears: TLabeledEdit
    Left = 379
    Top = 40
    Width = 50
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 58
    EditLabel.Height = 16
    EditLabel.Caption = 'Nr. Years:'
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
    TabOrder = 2
    Text = '50'
  end
  object edtRep: TLabeledEdit
    Left = 379
    Top = 10
    Width = 50
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 84
    EditLabel.Height = 16
    EditLabel.Caption = 'Nr. Replicates:'
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
    TabOrder = 1
    Text = '1'
    OnChange = edtRepChange
  end
  object edtSimNr: TLabeledEdit
    Left = 175
    Top = 10
    Width = 50
    Height = 24
    Alignment = taRightJustify
    EditLabel.Width = 152
    EditLabel.Height = 16
    EditLabel.Caption = 'Simulation number:       '
    EditLabel.Font.Charset = DEFAULT_CHARSET
    EditLabel.Font.Color = clWindowText
    EditLabel.Font.Height = -13
    EditLabel.Font.Name = 'Tahoma'
    EditLabel.Font.Style = [fsBold]
    EditLabel.ParentFont = False
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = []
    LabelPosition = lpLeft
    ParentFont = False
    TabOrder = 0
    Text = '0'
  end
  object BtnSeeding: TBitBtn
    Left = 451
    Top = 8
    Width = 209
    Height = 32
    Caption = 'Set Initialisation Rules'
    Font.Charset = DEFAULT_CHARSET
    Font.Color = clWindowText
    Font.Height = -13
    Font.Name = 'Tahoma'
    Font.Style = [fsBold]
    Glyph.Data = {
      36030000424D3603000000000000360000002800000010000000100000000100
      18000000000000030000120B0000120B00000000000000000000FF00FFFF00FF
      FF00FFFF00FFFF00FF044906055B09066C0C066C0C055E0A044C06FF00FFFF00
      FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF05600905600908911309B01809
      B31A09B31909B11907961405680C05680CFF00FFFF00FFFF00FFFF00FFFF00FF
      0A6A150A7F150BB61C09B91A08B41807B21609B31909B41909B81A09B91A0783
      10044D06FF00FFFF00FFFF00FF0B6A150F852216BD3411B7270BB21C07B11608
      B11709B21909B21909B21909B41909BA1A07841006670CFF00FFFF00FF0B6A15
      20BE491BBD4014B7300AB21F28BC36DFF5E1EEFAEF63CE6D09B21909B21909B3
      1909BA1A06670CFF00FF0872101B9A3A2AC65B1DBB450EB4250BB31B11B4219A
      DFA0FFFFFFF7FDF85ACB6509B21909B21909B81A089413045D090872102AB65B
      2CC56522BD4D0FB4220AB21A0CB31C0AB2198DDB95FDFEFDF6FCF758CB6309B2
      1909B51A08AB17045D090F821C37C26C33C76CCDF1DAC9EFD3C7EED0C8EFD2C5
      EED0C7EECFF8FDF9FFFFFFF2FBF36FD27908B41909B31905650B138D2358CC83
      42C977FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFDFEFDFFFFFFFFFFFFBCEA
      C10AB41A09B319066D0D0F911D6FD2935FD38D6DD49572D69971D69872D69964
      D28C92DFA8FBFEFBFFFFFFACE5B82EBF4C11B82B08B11905610A0F911D67CC83
      9BE5BA38C67030C36938C56F38C56F70D697E8F8EEFFFFFF9FE2B120BD481AB9
      3E10BA2908A31705610AFF00FF25AE39BCEDD282DBA428C0632FC26753CD82F7
      FDF9FFFFFF9CE2B222BC4B1DBA4118B73614C0300A8517FF00FFFF00FF25AE39
      71D28CD2F4E180DAA336C46D39C56FBCECCEABE6C22DC26324BE5623BC4D1FC1
      4616AE340A8517FF00FFFF00FFFF00FF25AE3984D89FDBF7EAAFE8C66BD49352
      CC8144C97849CA7B48CB7839CB6A21B6490F7C1FFF00FFFF00FFFF00FFFF00FF
      FF00FF25AE3925AE39ADE8C5CCF2DEBAEDD1A6E7C291E2B364D4922FB1572FB1
      57FF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FFFF00FF32B74E25AE3925
      AE3925AE3925AE3924A342FF00FFFF00FFFF00FFFF00FFFF00FF}
    ParentFont = False
    TabOrder = 3
    TabStop = False
    OnClick = BtnSeedingClick
  end
end
