define([
    'app',
    'app/model/thread/invite'
], function (App) {

    App.View.Thread.InviteCreateForm = App.Backbone.Form.extend({

        schema: {
            to_user: {
                type: 'Select',
                title: 'Who',
                options: function(callback) {
                    var users = new App.Collection.Users();

                    // TODO users.url
                    users
                        .fetch()
                        .then(function() {
                            var options = [];
                            users.each(function (user) {
                                var label = String(user.get('first_name') + ' ' + user.get('last_name')).trim();
                                if (0 == label.length) {
                                    label = '&lt;' + user.get('email') + '&gt;';
                                }
                                options.push({val: user.get('id'), label: label});
                            });
                            callback(options);
                        }, function(jqXHR, textStatus, errorThrown) {
                            throw textStatus + ' ' + errorThrown;
                        });
                },
                validators: ['required', App.Model.Invite.validators.to_user]
            }
        },

        template: App._.template('\
            <form class="thread-create-form form-horizontal" role="form">\
                <div data-fieldsets></div>\
                <div class="form-group form-footer buttons">\
                    <div class="col-sm-offset-2 col-sm-10">\
                        <button type="submit" class="btn btn-primary form-submit">\
                            Invite\
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

            var errors = this.commit();
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


    App.View.Thread.InviteCreateForm.create = function(options, createOptions) {
        options = options || {};
        if (!options.model) {
            options.model = App.Model.Invite.createNew();
        }

        var form = new App.View.Thread.InviteCreateForm(options);
        if (createOptions) {
            App.View.create(form, createOptions);
        }

        return form;
    };

    return App.View.Thread.InviteCreateForm;
});
