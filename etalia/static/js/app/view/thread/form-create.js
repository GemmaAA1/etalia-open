define([
    'underscore',
    'app',
    'app/collection/library/paper',
    'app/model/thread/post'
], function (_, App) {

    App.View.Thread.CreateForm = App.Backbone.Form.extend({

        schema: {
            title: {type: 'Text', validators: ['required']},

            paper: {type: 'Select', options: function(callback) {
                var user = App.Model.User.getCurrent(),
                    papers = new App.Collection.Papers();

                papers.url = App.config.api_root + '/user/user-libs/' + user.get('id') + '/papers';
                papers
                    .fetch()
                    .then(function() {
                        callback(papers);
                    }, function(jqXHR, textStatus, errorThrown) {
                        throw textStatus + ' ' + errorThrown;
                    });

            }, validators: ['required']}
        },

        template: _.template('\
            <form class="thread-create-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <div class="form-group form-footer">\
                    <div class="col-sm-offset-2 col-sm-10">\
                        <button type="submit" class="btn btn-primary form-submit">\
                            Create\
                        </button>\
                        <button type="button" class="btn btn-default form-cancel">\
                            Cancel\
                        </button>\
                    </div>\
                </div>\
            </form>\
        '),

        events: {
            "click .form-submit": "_onSubmit",
            "click .form-cancel": "_onCancel"
        },

        initialize: function () {
            App.Backbone.Form.prototype.initialize.apply(this, arguments);

            this.listenTo(this, "submit", this._onSubmit);
        },

        _onSubmit: function(e) {
            e.preventDefault();

            var errors = this.commit({ validate: true });
            if (errors) {
                App.log('Errors', errors);

                this.trigger('validation_failure', this);
                return;
            }

            this.trigger('validation_success', this);
        },

        _onCancel: function(e) {
            e.preventDefault();

            this.trigger('cancel', this);
        }
    });


    App.View.Thread.CreateForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = App.Model.Thread.createNew();
        }

        var form = new App.View.Thread.CreateForm(options);
        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.CreateForm;
});
