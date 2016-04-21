define([
    'require',
    'underscore',
    'backbone',
    'handlebars',
    'moment',
    'app/ui/layout',
    'app/ui/controls',
    'backbone/relational',
    'backbone/paginator',
    'backbone/forms'
], function(require, _, Backbone, Handlebars, Moment, Layout, Controls) {

    /**
     * Backbone
     */
    var _sync = Backbone.sync;
    Backbone.sync = function (method, model, options) {
        // Add trailing slash to backbone model views
        var _url = _.isFunction(model.url) ? model.url() : model.url;
        _url += _url.charAt(_url.length - 1) == '/' ? '' : '/';
        options = _.extend(options, {
            url: _url
        });
        return _sync(method, model, options);
    };
    Backbone.Collection.prototype.parse = function (response, options) {
        if (_.has(response, 'results')) {
            return response.results;
        }
        return response;
    };


    /**
     * Backbone Relational
     */
    function isApiUrl(url) {
        return _.isString(url) && /http:\/\/[0-9a-z-:\.]+\/[a-z0-9-\/]+\/[0-9]+\//.test(url);
    }

    function extractIdFromApiUrl(url) {
        // TODO assert that the ID is properly extracted.
        return parseInt(/\d+/.exec(/\/\d+\//.exec(url)[0])[0]);
    }

    var relationalModel_updateRelations =
        Backbone.RelationalModel.prototype.__updateRelations =
            Backbone.RelationalModel.prototype.updateRelations;
    Backbone.RelationalModel.prototype.updateRelations = function (changedAttrs, options) {
        if (this._isInitialized && !this.isLocked()) {
            _.each(this._relations, function (rel) {
                var value;
                if (rel.keySource in changedAttrs) {
                    value = this.attributes[rel.keySource];
                    if (isApiUrl(value)) {
                        this.attributes[rel.keySource] = extractIdFromApiUrl(value);
                        changedAttrs[rel.keySource] = extractIdFromApiUrl(value);
                    }
                }
                if (rel.key in changedAttrs) {
                    value = this.attributes[rel.key];
                    if (isApiUrl(value)) {
                        this.attributes[rel.key] = extractIdFromApiUrl(value);
                        changedAttrs[rel.key] = extractIdFromApiUrl(value);
                    }
                }
            }, this);
        }

        relationalModel_updateRelations.call(this, changedAttrs, options);
    };

    var hasOneSetKeyContents = Backbone.HasOne.prototype.__setKeyContents = Backbone.HasOne.prototype.setKeyContents;
    Backbone.HasOne.prototype.setKeyContents = function( keyContents ) {
        // Extract ID from Api Url
        if (isApiUrl(keyContents)) {
            keyContents = extractIdFromApiUrl(keyContents);
        }
        hasOneSetKeyContents.call( this, keyContents );
    };

    // @TODO test
    var hasManySetKeyContents = Backbone.HasMany.prototype.__setKeyContents = Backbone.HasMany.prototype.setKeyContents;
    Backbone.HasMany.prototype.setKeyContents = function( keyContents ) {

        // Asserts that keyContents is an array of api url string
        if (_.isArray(keyContents)) {
            var assert = true;
            _.each(keyContents, function(value) {
                if (!isApiUrl(value)) {
                    assert = false;
                }
            });
        }

        // Convert url[] to id[]
        if (assert) {
            keyContents = _.map(keyContents, function(value) {
                return extractIdFromApiUrl(value);
            });
        }

        hasManySetKeyContents.call( this, keyContents );
    };


    /**
     * Backbone Paginator
     */
    Backbone.PageableCollection.prototype.parseRecords = function (response) {
        if (_.has(response, 'results')) {
            return response.results;
        }
        return response;
    };
    Backbone.PageableCollection.prototype.parseLinks = function (response) {
        var links = {};
        if (_.has(response, 'next')) {
            links.next = response.next;
        }
        if (_.has(response, 'prev')) {
            links.prev = response.prev;
        }
        return links;
    };


    /**
     * Backbone Forms
     */
        // @see https://github.com/gerardobort/bbf-bbr-integration/blob/master/test.js
    Backbone.Form.editors.Base.prototype._initialize = Backbone.Form.editors.Base.prototype.initialize;
    Backbone.Form.editors.Base.prototype.initialize = function (options) {
        Backbone.Form.editors.Base.prototype._initialize.call(this, options);
        // this patch adds compatibility between backbone-forms and backbone-relational libraries
        if (options.model instanceof Backbone.RelationalModel && options.model.get(options.key) instanceof Backbone.Collection) {
            this.value = options.model.get(options.key).toJSON();
        }
    };

    Backbone.Form.template = _.template('\
    <form class="form-horizontal" role="form">\
      <div data-fieldsets></div>\
      <% if (submitButton) { %>\
        <button type="submit" class="btn btn-primary"><%= submitButton %></button>\
      <% } %>\
    </form>\
  ');
    Backbone.Form.Field.template = _.template('\
    <div class="form-group field-<%= key %>">\
      <% if (title) { %>\
      <label class="col-sm-2 control-label hidden-xs" for="<%= editorId %>"><%= title %></label>\
      <% } %>\
      <div class="col-sm-<% if (title) { %>10<% } else { %>12<% } %>">\
        <span data-editor></span>\
        <p class="help-block" data-error></p>\
        <p class="help-block"><%= help %></p>\
      </div>\
    </div>\
  ');

    Backbone.Form.editors.Radio = Backbone.Form.editors.Radio.extend({
        tagName: 'div',
        className: 'btn-group',
        attributes: {
            'data-toggle': 'buttons'
        }
    }, {
        template: _.template('\
    <% _.each(items, function(item) { %>\
      <label for="<%= item.id %>" class="btn btn-default">\
        <input type="radio" name="<%= item.name %>" value="<%= item.value %>" id="<%= item.id %>" autocomplete="off"> <%= item.label %>\
      </label>\
    <% }); %>\
  ', null, Backbone.Form.templateSettings)
    });

    Backbone.Form.editors.Tinymce = Backbone.Form.editors.TextArea.extend({

        /*tagName: 'textarea',

        events: {
            'change': function() {
                // The 'change' event should be triggered whenever something happens
                // that affects the result of `this.getValue()`.
                this.trigger('change', this);
            },
            'focus': function() {
                // The 'focus' event should be triggered whenever an input within
                // this editor becomes the `document.activeElement`.
                this.trigger('focus', this);
                // This call automatically sets `this.hasFocus` to `true`.
            },
            'blur': function() {
                // The 'blur' event should be triggered whenever an input within
                // this editor stops being the `document.activeElement`.
                this.trigger('blur', this);
                // This call automatically sets `this.hasFocus` to `false`.
            }
        },*/

        defaultSettings: {
            resize: false,
            schema: "html5-strict",
            fix_list_elements : true,
            keep_styles: false,
            invalid_elements : "span",
            statusbar: false,
            menubar: false,
            height: 160,
            plugins: "lists,paste,link,searchreplace,hr",
            toolbar: "undo redo | formatselect | bold italic | bullist numlist outdent indent | link hr blockquote",
            skin_url: '/static/css/lib/tinymce',
            body_class: 'element-content',
            content_css: [
                'https://fonts.googleapis.com/css?family=Lato:400,700,700italic,400italic',
                '/static/css/lib/bootstrap.css',
                '/static/css/app/content.css'
            ]
        },

        tinymceEditorId: null,

        initialize: function(options) {
            Backbone.Form.editors.Base.prototype.initialize.call(this, options);
        },

        render: function() {
            this.setValue(this.value);

            var that = this,
                settings = {};

            this.tinymceEditorId = that.$el.attr('id');

            _.defaults(
                settings,
                this.schema.settings,
                this.defaultSettings,
                {selector: '#' + that.tinymceEditorId}
            );
            settings.setup = function(editor) {
                var onChange = function() {
                    that.setValue(editor.getContent());

                    /*und.convert(editor.getContent(), function(err, markdown) {
                     $('#markdown').val(markdown);
                     }, { keepHtml: false });*/
                };

                editor.on('keyup', onChange);
                editor.on('change', onChange);
            };

            require(['jquery', 'tinymce'], function($, tinymce) {
                tinymce.EditorManager.execCommand('mceRemoveEditor', true, that.tinymceEditorId);

                tinymce.init(settings);
            });

            return this;
        },

        remove: function() {
            var that = this;
            require(['tinymce'], function(tinymce) {
                tinymce.EditorManager.execCommand('mceRemoveEditor', true, that.tinymceEditorId);
            });

            if (this.tinymceEditor) {
                this.tinymceEditor.destroy();
                this.tinymceEditor = null;
            }

            Backbone.View.prototype.remove.call(this);
        }

        /*getValue: function() {
            return this.$el.val();
        },

        setValue: function(value) {
            this.$el.val(value);
        },

        focus: function() {
            if (this.hasFocus) return;

            // This method call should result in an input within this editor
            // becoming the `document.activeElement`.
            // This, in turn, should result in this editor's `focus` event
            // being triggered, setting `this.hasFocus` to `true`.
            // See above for more detail.
            this.$el.focus();
        },

        blur: function() {
            if (!this.hasFocus) return;

            this.$el.blur();
        }*/
    });


    /**
     * Handlebars helpers
     */
    Handlebars.registerHelper('ifCond', function (a, operator, b, options) {
        switch (operator) {
            case '==':
                return (a == b) ? options.fn(this) : options.inverse(this);
            case '===':
                return (a === b) ? options.fn(this) : options.inverse(this);
            case '<':
                return (a < b) ? options.fn(this) : options.inverse(this);
            case '<=':
                return (a <= b) ? options.fn(this) : options.inverse(this);
            case '>':
                return (a > b) ? options.fn(this) : options.inverse(this);
            case '>=':
                return (a >= b) ? options.fn(this) : options.inverse(this);
            case 'and':
                return (a && b) ? options.fn(this) : options.inverse(this);
            case 'andNot':
                return (a && !b) ? options.fn(this) : options.inverse(this);
            case 'or':
                return (a || b) ? options.fn(this) : options.inverse(this);
            case 'orNot':
                return (a || !b) ? options.fn(this) : options.inverse(this);
            default:
                return options.inverse(this);
        }
    });
    Handlebars.registerHelper('thread_pin_class', function() {
        if (this.state && this.state.get('watch') === App.Model.State.WATCH_PINNED) {
            return ' active';
        }
        return '';
    });
    Handlebars.registerHelper('thread_ban_class', function() {
        if (this.state && this.state.get('watch') === App.Model.State.WATCH_BANNED) {
            return ' active';
        }
        return '';
    });
    Handlebars.registerHelper('full_name', function(user) {
        /*if (!user) {
            return '';
        }*/
        return user.get('first_name') + " " + user.get('last_name');
    });
    Handlebars.registerHelper('paper_title_authors', function(paper) {
        /*if (!paper) {
            return '';
        }*/
        var authors = paper.get('authors').map(function (author) {
            return author.get('first_name') + " " + author.get('last_name');
        });
        var output = paper.get('title') + " (" + authors.splice(0, 4).join(', ') + ")";

        return new Handlebars.SafeString(output);
    });
    Handlebars.registerHelper('date', function(date) {
        return Moment(date).format('MMM D, YYYY');
    });


    /**
     * App
     */
    window.App = App = {
        Backbone: Backbone,
        Handlebars: Handlebars,

        Layout: Layout,
        Controls: Controls,

        Const: {},
        Model: {},
        Collection: {},
        View: {
            createDefaults: {
                $target: null,
                append: false
            },
            create: function(view, options) {
                options = _.defaults(options, App.View.createDefaults);
                if (options.$target && options.$target.length) {
                    if (true === options.append) {
                        options.$target.append(view.render().$el);
                    } else {
                        options.$target.replaceWith(view.render().$el);
                    }
                }
            }
        },

        config: {
            debug: true,
            api_root: '/api/v1'
        },

        defaults: function(options, defaults) {
            return _.extend(options, defaults);
        },

        log: function() {
            if (this.config.debug) {
                if (arguments.length > 1) {
                    console.log(arguments[0], Array.prototype.splice.call(arguments, 1));
                } else {
                    console.log(arguments[0]);
                }
            }
        }
    };

    return App;
});
